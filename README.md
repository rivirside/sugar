# SUGAR v2

Systematic Utilization of Glycans for Alternate Routes -- a computational platform for discovering synthesis pathways between sugars.

## Structure

- `pipeline/` -- Python data pipeline (generates compound/reaction JSON)
- `web/` -- Next.js frontend (deployed to Vercel)

## Quick Start

### Generate data
```bash
cd pipeline/..
pip install -r pipeline/requirements.txt
python -m pipeline.run_pipeline
```

### Run frontend
```bash
cd web
npm install
npm run dev
```
