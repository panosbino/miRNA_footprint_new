#!/bin/bash
# find_and_restore.sh
# Searches git history for deleted files and offers to restore them

FILES=(
"resources/Targets_Negative_Control_mESCs.rds"
"resources/Targets__combined_124-3p_124-3p.2_506-3p.rds"
"resources/Targets__combined_124-3p_124-3p.2_506-3p_top100.rds"
"resources/Targets__combined_124-3p_124-3p.2_506-3p_top1000.rds"
"resources/Targets__combined_124-3p_124-3p.2_506-3p_top500.rds"
"resources/Targets__miR-124-3p.rds"
"resources/Targets__miR-124-3p_top100.rds"
"resources/Targets__miR-124-3p_top1000.rds"
"resources/Targets__miR-124-3p_top500.rds"
"resources/Targets_combined_miR_199-5p_miR_199-3p.rds"
"resources/Targets_combined_top3_families_mESCs.rds"
"resources/Targetscan_files/TargetScan8.0__miR-124-3p.1.Human.predicted_targets.txt"
"resources/Targetscan_files/TargetScan8.0__miR-124-3p.2_506-3p.Human.predicted_targets.txt"
"resources/Targetscan_files/human/TargetScan8.0__miR-199-3p.predicted_targets.txt"
"resources/Targetscan_files/human/TargetScan8.0__miR-199-5p.predicted_targets.txt"
"resources/Targetscan_files/mouse/TargetScan8.0__miR-1-3p_206-3p.predicted_targets.txt"
"resources/Targetscan_files/mouse/TargetScan8.0__miR-122-5p.predicted_targets.txt"
"resources/Targetscan_files/mouse/TargetScan8.0__miR-17-5p_20-5p_93-5p_106-5p.predicted_targets.txt"
"resources/Targetscan_files/mouse/TargetScan8.0__miR-291-3p_294-3p_295-3p_302-3p.predicted_targets.txt"
"resources/Targetscan_files/mouse/TargetScan8.0__miR-292a-3p_467a-5p.predicted_targets.txt"
"resources/Targetscan_files/mouse/TargetScan8.0__miR-9-5p.predicted_targets.txt"
)

echo "=== Step 1: Finding the last commit where each file existed ==="
echo ""

RESTORE_LOG="restore_commands.sh"
echo "#!/bin/bash" > "$RESTORE_LOG"
echo "# Auto-generated restore commands - review before running" >> "$RESTORE_LOG"
echo "" >> "$RESTORE_LOG"

for f in "${FILES[@]}"; do
    echo "--- $f ---"
    # Find the last commit that touched this file (across all history, even unreachable)
    commit=$(git log --all --diff-filter=D --pretty=format:"%H" -- "$f" | head -n 1)

    if [ -z "$commit" ]; then
        # fallback: check dangling objects too
        commit=$(git rev-list --objects --all 2>/dev/null | grep -F "$f" | awk '{print $1}' | head -n 1)
    fi

    if [ -n "$commit" ]; then
        # get the parent commit (state BEFORE deletion)
        parent="${commit}^"
        echo "  Found deletion commit: $commit"
        echo "  Restoring from parent: $parent"
        mkdir_cmd="mkdir -p \"$(dirname "$f")\""
        restore_cmd="git show ${parent}:\"$f\" > \"$f\""
        echo "$mkdir_cmd" >> "$RESTORE_LOG"
        echo "$restore_cmd" >> "$RESTORE_LOG"
    else
        echo "  NOT FOUND in history — may need to check Dardel or backups"
    fi
    echo ""
done

chmod +x "$RESTORE_LOG"
echo "=== Done. Review '$RESTORE_LOG' before running it. ==="
echo "Run: ./$RESTORE_LOG"
