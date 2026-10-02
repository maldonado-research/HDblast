#!/usr/bin/env python3
"""Read-only AST proof: only NPZ storage/lifetime/hash implementation changes."""
from __future__ import annotations
import ast
import copy
import hashlib
import json
from pathlib import Path


def require(condition,message):
    if not condition:
        raise RuntimeError(message)


def dump(node):
    return ast.dump(node,include_attributes=False)


def called(node,owner,name):
    return (isinstance(node,ast.Call) and isinstance(node.func,ast.Attribute) and
            isinstance(node.func.value,ast.Name) and node.func.value.id==owner and node.func.attr==name)


class StorageOnly(ast.NodeTransformer):
    def __init__(self):
        self.removed=[]

    def visit_With(self,node):
        if (len(node.items)==1 and isinstance(node.items[0].context_expr,ast.Call) and
            isinstance(node.items[0].context_expr.func,ast.Name) and node.items[0].context_expr.func.id=='StreamingNpz' and
            isinstance(node.items[0].optional_vars,ast.Name) and node.items[0].optional_vars.id=='archive'):
            self.removed.append('StreamingNpz archive lifetime context')
            result=[]
            for item in node.body:
                visited=self.visit(item)
                if visited is not None:
                    result.extend(visited if isinstance(visited,list) else [visited])
            return result
        return self.generic_visit(node)

    def visit_Assign(self,node):
        if len(node.targets)==1 and isinstance(node.targets[0],ast.Name):
            if node.targets[0].id=='archive' and isinstance(node.value,ast.Dict):
                require([key.value for key in node.value.keys]==['k','momentum_weights'],'Unexpected archive grid initialization')
                self.removed.append('in-memory archive grid dict')
                return None
            if node.targets[0].id=='path':
                self.removed.append('NPZ output path placement')
                return None
        return self.generic_visit(node)

    def visit_Expr(self,node):
        call=node.value
        if called(call,'archive','close'):
            require(not call.args and not call.keywords,'Unexpected archive close arguments')
            self.removed.append('stream close replaces bulk savez')
            return None
        if called(call,'np','savez_compressed'):
            require(dump(call)==dump(ast.parse('np.savez_compressed(path,**archive)').body[0].value),
                    'Unexpected prior bulk archive call')
            self.removed.append('bulk savez-compressed archive serialization')
            return None
        if called(call,'archive','update') and len(call.args)==1 and isinstance(call.args[0],ast.Dict):
            keys=call.args[0].keys
            if all(isinstance(key,ast.Constant) for key in keys) and [key.value for key in keys]==['k','momentum_weights']:
                self.removed.append('streamed archive grid initialization')
                return None
        if called(call,'active_arrays','update') and len(call.args)==1 and isinstance(call.args[0],ast.DictComp):
            key=call.args[0].key
            require(isinstance(key,ast.BinOp) and isinstance(key.left,ast.Constant) and key.left.value=='archived_',
                    'Unexpected active-array mutation')
            self.removed.append('retained archived_* ndarray aliases')
            return None
        return self.generic_visit(node)

    def visit_Delete(self,node):
        if [target.id for target in node.targets if isinstance(target,ast.Name)]==['contact_archive','integrands']:
            self.removed.append('release streamed observation-only dictionaries')
            return None
        return self.generic_visit(node)


def main():
    here=Path(__file__).resolve().parent
    prior_path=here/'PRIOR_FAILED_PRODUCER.py.txt'
    current_path=here/'forced_metric.py'
    prior_text,current_text=prior_path.read_text(),current_path.read_text()
    prior,current=ast.parse(prior_text),ast.parse(current_text)
    before={node.name:node for node in prior.body if isinstance(node,ast.FunctionDef)}
    after={node.name:node for node in current.body if isinstance(node,ast.FunctionDef)}
    require(set(before)==set(after),'Producer function inventory changed')
    identical=[]
    for name in before:
        if name not in ('sha','provenance_gate','evolve'):
            require(ast.get_source_segment(prior_text,before[name])==ast.get_source_segment(current_text,after[name]),
                    'Non-archive producer function bytes changed: '+name)
            identical.append(name)
    constants_before=[dump(node) for node in prior.body if isinstance(node,ast.Assign)]
    constants_after=[dump(node) for node in current.body if isinstance(node,ast.Assign)]
    require(constants_before==constants_after,'Model/grid/gate/budget declarations changed')
    old,new=StorageOnly(),StorageOnly()
    old_ast=old.visit(copy.deepcopy(before['evolve']))
    new_ast=new.visit(copy.deepcopy(after['evolve']))
    require(dump(old_ast)==dump(new_ast),'A numerical/evolution/archive-field expression changed')
    old_gate,new_gate=copy.deepcopy(before['provenance_gate']),copy.deepcopy(after['provenance_gate'])
    additions=0
    for node in ast.walk(new_gate):
        if isinstance(node,ast.Set):
            kept=[item for item in node.elts if not (isinstance(item,ast.Constant) and item.value=='stream_npz.py')]
            additions+=len(node.elts)-len(kept);node.elts=kept
    require(additions==1 and dump(old_gate)==dump(new_gate),'Provenance changes exceed the mandatory stream module pin')
    require(dump(after['sha'].body[0])==dump(ast.parse('return file_sha256(path)').body[0]),'Unexpected bounded SHA wrapper')
    sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
    report={'schema_version':1,'status':'PASS','scope':'Read-only source/AST identity; no source values, modes, stresses or physical integrals evaluated.',
            'physical_evaluations':0,'prior_failed_producer_sha256':sha(prior_path),'new_producer_sha256':sha(current_path),
            'byte_identical_functions':identical,'model_grid_gates_budgets_exact_AST_identity':True,
            'evolution_numerical_and_retained_field_AST_identity':True,'old_storage_nodes_removed':old.removed,
            'new_storage_nodes_removed':new.removed,'provenance_only_addition':'stream_npz.py must be present in frozen manifest',
            'hash_only_change':'bounded 1MiB SHA256 chunks; identical digest',
            'stream_module_sha256':sha(here/'stream_npz.py'),'audit_source_sha256':sha(Path(__file__))}
    (here/'ARCHIVE_ONLY_INVARIANCE.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
