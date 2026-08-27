from __future__ import annotations
import argparse, json
from .real_intervention_bundle import ingest_measured_study

def main():
 p=argparse.ArgumentParser(description='Strict measured intervention study ingestion')
 p.add_argument('--data',required=True);p.add_argument('--manifest',required=True);p.add_argument('--out',required=True);p.add_argument('--require-files',action='store_true');a=p.parse_args()
 r=ingest_measured_study(a.data,a.manifest,a.out,require_files=a.require_files)
 print(json.dumps(r['ingest'],indent=2)); return 0 if r['ingest']['status']=='accepted_for_locked_analysis' else 2
if __name__=='__main__': raise SystemExit(main())
