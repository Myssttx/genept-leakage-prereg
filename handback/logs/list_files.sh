#!/bin/bash
# Lists new/changed non-raw-data files with sizes (bytes) and SHA-256, for handback section 8.
cd ~/Documents/genept-leakage
echo '| Path | Size (bytes) | SHA-256 |'
echo '|---|---|---|'
for f in CLAUDE.md STAGES.md review.md plan.md requirements.txt .gitignore code/stage1_data_report.py code/tests/test_stage1_report.py data/README.md \
         handback/stage1_data_report.txt handback/stage1_pip_freeze.txt handback/stage1_pytest.txt handback/stage1_step6_models.txt \
         $(ls handback/logs/* | sort) $(ls data/papers/* | sort) $(ls data/web/*.* data/web/ncbi/* | sort); do
  [ -f "$f" ] && printf '| `%s` | %s | `%s` |\n' "$f" "$(stat -f %z "$f")" "$(shasum -a 256 "$f" | cut -d' ' -f1)"
done
echo
echo "Directories: \`handback/stage1_code_before/\` (copy of \`code/\` before the stage), \`data/models/\` ($(du -sk data/models | cut -f1) KiB, Hugging Face cache), \`.venv/\` ($(du -sk .venv | cut -f1) KiB)."
