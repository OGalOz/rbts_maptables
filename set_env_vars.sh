#!/bin/bash


chome=$PWD
impl=$PWD/lib/rbts_maptables/rbts_maptablesImpl.py
tst=$PWD/test/rbts_maptables_server_test.py
util_dir=$PWD/lib/util
tmp_dir=$PWD/test_local/workdir/tmp
ui_dir=$PWD/ui/narrative/methods/run_rbts_maptables/
uspec=$ui_dir/spec.json
udisp=$ui_dir/display.yaml
Trash=$PWD/../Trash

#Docker fix
docker run -it -v /var/run/docker.sock:/run/docker.sock alpine chmod g+w /run/docker.sock

# clean up
find . -name '.DS_Store' -type f -delete
