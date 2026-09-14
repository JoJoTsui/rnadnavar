#!/usr/bin/env python3
"""Create a content-identity manifest for frozen policy validation inputs."""
import argparse, hashlib, json
from pathlib import Path

def digest(path):
 h=hashlib.sha256();
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
 return h.hexdigest()

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--out',type=Path,required=True); ap.add_argument('paths',nargs='+',type=Path); a=ap.parse_args()
 rows=[]; missing=[]
 for path in a.paths:
  path=path.resolve()
  if not path.is_file(): missing.append(str(path)); continue
  rows.append({"path":str(path),"sha256":digest(path),"size_bytes":path.stat().st_size})
 result={"status":"pass" if not missing else "blocked","immutable":True,"files":rows,"missing":missing}
 a.out.parent.mkdir(parents=True,exist_ok=True); a.out.write_text(json.dumps(result,indent=2)+"\n"); print(json.dumps({"status":result["status"],"files":len(rows),"missing":len(missing)})); return 0 if not missing else 1
if __name__=='__main__': raise SystemExit(main())
