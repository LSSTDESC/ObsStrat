Building the paper requires checking out the des-tex submodule.

So, the complete sequence to clone the git repo and build the paper is
as follows:

```sh
git clone git@github.com:LSSTDESC/obs_strat.git obs_strat
cd obs_strat
git checkout issue/12/writeup
git submodule init
git submodule update
cd writeup
make all
```

`ditaa` was used to create the wedding-cake diagram. It can be
regenerated thus:

```
java -jar ~/ditaa/ditaa0_9.jar -E -S wedding_cake.ditaa wedding_cake.png
```
