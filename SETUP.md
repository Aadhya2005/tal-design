# Running this on your own machine

## 1. Check Python

```bash
python3 --version
```

Needs 3.10 or newer. If that command is not found, install Python from
python.org (tick "Add Python to PATH" on Windows) and reopen the terminal.

## 2. Open it in VS Code (optional, but recommended)

Install [VS Code](https://code.visualstudio.com) and the "Python" extension
(Extensions icon in the sidebar -> search "Python" -> install the one by
Microsoft). Then File -> Open Folder -> this project's folder.

The sidebar now shows the project as a file tree, and `.vscode/settings.json`
already points VS Code at the virtual environment you're about to create
below, so once it exists VS Code picks it up automatically. Use
Terminal -> New Terminal for the commands in the steps below - that opens a
terminal panel inside the editor rather than a separate window.

## 3. Set up a virtual environment

From inside the project folder (in VS Code's terminal panel, or any terminal):

```bash
python3 -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

A virtual environment keeps these packages out of your system Python. Your
prompt should now start with `(.venv)`. If VS Code asks you to select an
interpreter, or shows a notification about one, pick the one inside `.venv`.

## 4. Put the input files in place

Neither input is committed to the repo. Download them yourself:

| File | Where | Save as |
|---|---|---|
| RsTAL protein | uniprot.org/uniprotkb/Q3IWB0 -> Download -> Text format | `data/Q3IWB0.txt` |
| pET-28a(+) | SnapGene plasmid library or Addgene | `data/pET-28a.dna` or `.gb` |

## 5. Run it

```bash
python3 scripts/run_all.py --protein data/Q3IWB0.txt --vector data/pET-28a.dna --seed 1
python3 scripts/audit.py
```

Expect 8/8 checks from the first command and ALL CHECKS PASSED from the second.
Everything lands in `output/`.

In VS Code you can instead open the Run and Debug panel (the triangle-with-bug
icon in the left sidebar) and run the two pre-configured entries, "Run
pipeline (run_all.py)" and "Run independent audit (audit.py)", which does the
same thing without typing the commands.

Once it runs, open `output/` in the sidebar: the `.png` figures preview
directly in VS Code, and right-clicking `DESIGN_REPORT.md` -> Open Preview
renders it instead of showing raw Markdown.

## 6. If something breaks

| Symptom | Cause | Fix |
|---|---|---|
| `ModuleNotFoundError` | venv not active | re-run the `activate` line |
| `T7 promoter not found on either strand` | wrong vector file | check you downloaded pET-28a(+), not another pET |
| `expected exactly 1 NdeI site, found 0` | vector already has an insert | use the empty vector |
| Different numbers from the report | seed not passed | add `--seed 1` |
| Optimisation takes several minutes | normal | the solver is doing real work |

## 7. Change something

The point of having the code is that you can. A few one-line experiments:

- `--seed 7` gives a different valid optimised sequence. Compare the CAI.
- In `scripts/optimise.py`, widen `ramp_gc` to `(0.25, 0.70)` and re-run. Watch
  the forward primer's hairpin dG get worse. That is the 5' ramp earning its place.
- Add an enzyme to `FORBIDDEN_SITES` and confirm it disappears from the insert.

Record what you tried and what happened. That is the part of a project that is
actually yours.
