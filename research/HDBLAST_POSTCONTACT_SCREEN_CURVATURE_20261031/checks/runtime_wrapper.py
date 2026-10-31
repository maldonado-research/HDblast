#!/usr/bin/env python3
"""Override runtime path globals without changing retained scientific sources."""
import argparse
from pathlib import Path
import runpy
import sys

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--script',required=True,type=Path)
    p.add_argument('--set-path',nargs=2,action='append',default=[])
    p.add_argument('script_args',nargs=argparse.REMAINDER)
    a=p.parse_args()
    namespace=runpy.run_path(str(a.script),run_name='portable_control_module')
    entry=namespace.get('main')
    if not callable(entry):raise RuntimeError('RUNTIME_MAIN_REQUIRED')
    for key,value in a.set_path:
        if key not in ('SOURCE','INPUT'):raise RuntimeError('RUNTIME_PATH_OVERRIDE_ONLY')
        if key not in entry.__globals__:raise RuntimeError('RUNTIME_PATH_GLOBAL_MISSING')
        entry.__globals__[key]=Path(value)
    args=a.script_args[1:] if a.script_args[:1]==['--'] else a.script_args
    sys.argv=[str(a.script)]+args
    entry()

if __name__=='__main__':main()
