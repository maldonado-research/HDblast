from pathlib import Path
import argparse,os,runpy,sys
parser=argparse.ArgumentParser();parser.add_argument('--proof',required=True);parser.add_argument('--output',required=True);args=parser.parse_args()
proof=Path(args.proof).resolve();root=proof.parent
sys.dont_write_bytecode=True
sys.path.insert(0,str(root))
write_flags=os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC|os.O_APPEND
def audit(event,values):
 if event!='open':return
 file,mode,flags=values
 if not isinstance(file,(str,bytes,os.PathLike)):return
 writing=(isinstance(mode,str) and any(c in mode for c in 'wax+')) or (isinstance(flags,int) and flags&write_flags)
 if writing and Path(file).resolve().is_relative_to(root):
  raise RuntimeError('Proof attempted to write its source tree: '+str(file))
sys.addaudithook(audit)
sys.argv=[str(proof),'--output',args.output]
runpy.run_path(str(proof),run_name='__main__')
