from __future__ import annotations
import argparse, json
from .flagship_benchmark import write_report

def main():
 p=argparse.ArgumentParser(description='Locked donor/clone/timepoint benchmark')
 p.add_argument('--input',required=True);p.add_argument('--out',required=True);a=p.parse_args()
 r=write_report(a.input,a.out); print(json.dumps({'out':a.out,'status':r['status'],'fingerprint':r['fingerprint']},indent=2))
if __name__=='__main__': main()
