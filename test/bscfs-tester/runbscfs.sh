FLOOR=/usr/local
$FLOOR/bb/scripts/fakelsf.sh \
--hostlist pike16 \
--ssd_min 200G \
--bscfs \
--pfs_path /mnt/scratch/pfs \
-- \
mpirun \
--allow-run-as-root \
--host pike16 --np 1 \
--x PFS_PATH --x BSCFS_PATH \
$FLOOR/bscfs/tests/chkpnt_write \
--compute_time \
--chkpnt_count 2 \
--header_size 4KiB \
--chkpnt_size 100GiB \
--stripe_size 512MiB \
--keep_all \
--chkpnt_dir bscfs_run
