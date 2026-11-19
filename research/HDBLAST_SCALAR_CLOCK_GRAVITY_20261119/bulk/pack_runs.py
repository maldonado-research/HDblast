"""Losslessly pack native arrays with byte shuffling and XORs against each matching zero run."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import zipfile
import numpy as np

HERE=Path(__file__).resolve().parent


def read_run(path):
    with np.load(path,allow_pickle=False) as data:
        return {key:data[key] for key in data}


def base_name(name):
    return "_".join(name.split("_")[:2])+"_zero"


def pack(source,archive):
    manifest={}
    zero={}
    def write(z,name,raw):
        info=zipfile.ZipInfo(name,(2026,11,19,0,0,0))
        info.compress_type=zipfile.ZIP_DEFLATED
        info.external_attr=0o100644<<16
        z.writestr(info,raw,compresslevel=7)
    names=sorted((p.stem for p in source.glob("*.npz")),key=lambda n:(not n.endswith("_zero"),n))
    with zipfile.ZipFile(archive,"w") as z:
        for name in names:
            data=read_run(source/(name+".npz"))
            base=zero.get(base_name(name),{})
            streams=[]
            layout=[]
            for key in sorted(data):
                a=data[key]
                raw=np.frombuffer(a.tobytes(),dtype=np.uint8)
                same=key in base and a.shape==base[key].shape and a.dtype==base[key].dtype
                if same:
                    raw=np.bitwise_xor(raw,np.frombuffer(base[key].tobytes(),dtype=np.uint8))
                streams.append(raw.reshape(-1,a.itemsize).T.tobytes())
                layout.append([key,a.dtype.str,list(a.shape),same])
            if name.endswith("_zero"):
                zero[name]=data
            raw=b"".join(streams)
            write(z,name+".bin",raw)
            manifest[name]={"layout":layout,"sha256":hashlib.sha256(raw).hexdigest()}
        write(z,"manifest.json",json.dumps(manifest,separators=(",",":")).encode())


def unpack(archive,target,check=None):
    target.mkdir(exist_ok=True,parents=True)
    with zipfile.ZipFile(archive) as z:
        manifest=json.loads(z.read("manifest.json"))
        zero={}
        for name,meta in sorted(manifest.items(),key=lambda row:(not row[0].endswith("_zero"),row[0])):
            raw=z.read(name+".bin")
            if hashlib.sha256(raw).hexdigest()!=meta["sha256"]:
                raise ValueError(name)
            base=zero.get(base_name(name),{})
            ix=0
            result={}
            for key,dtype,shape,same in meta["layout"]:
                dt=np.dtype(dtype)
                count=int(np.prod(shape,dtype=np.int64))*dt.itemsize
                shuffled=np.frombuffer(raw[ix:ix+count],dtype=np.uint8)
                arr=shuffled.reshape(dt.itemsize,-1).T.reshape(-1)
                if same:
                    arr=np.bitwise_xor(arr,np.frombuffer(base[key].tobytes(),dtype=np.uint8))
                result[key]=np.frombuffer(arr.tobytes(),dtype=dt).reshape(shape)
                ix+=count
            if ix!=len(raw):
                raise ValueError("Packed length mismatch "+name)
            if name.endswith("_zero"):
                zero[name]=result
            if check is not None:
                expected=read_run(check/(name+".npz"))
                if expected.keys()!=result.keys() or any(expected[k].tobytes()!=v.tobytes() for k,v in result.items()):
                    raise ValueError("Round trip changed native array "+name)
            else:
                np.savez_compressed(target/(name+".npz"),**result)
    return len(manifest)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("action",choices=("pack","unpack","verify"))
    ap.add_argument("--runs",type=Path,default=HERE/"runs")
    ap.add_argument("--archive",type=Path,default=HERE/"RUNS.zip")
    args=ap.parse_args()
    if args.action=="pack":
        pack(args.runs,args.archive)
    elif args.action=="unpack":
        print("Decoded runs",unpack(args.archive,args.runs))
    else:
        print("Bitwise round-tripped runs",unpack(args.archive,args.runs,check=args.runs))
