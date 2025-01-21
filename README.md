# tree_maker

This repository contains a macro which reads dst files and creates output root files containing proton candidate tracks relevant for the Azimuthal Partitioning analysis.

Once flow calculations were added to this analysis, the running was split into three parts:
1. Phi Flattening: This step produces histograms of the uncorrected phi distributions for each particle species. These histograms are then used to calculate the flow coefficients to flatten the phi distributions.
2. Event Plan Flattening: This step produces histograms of the uncorrected event plane distributions. Particle phi distributions are corrected with coefficients from the previous step. These histograms are then used to calculate the flow coefficients to flatten the event plane.
3. Tree Making: This step produces the final output trees of proton (pion) candidates. In addition, corrected event plane angles are included in the tree for each event.

### How to Run

To submit all jobs for a particular energy to the condor queue, edit the appropriate .xml file in the subs/ directory.
These submission xml files are very similar for all three steps of the analysis, so the final Tree Making step is used as an example.

What to change:

- `<!ENTITY energy "7">` - change to the desired energy
- `<!ENTITY prod "P10ih">` - change to the desired production
- `<!ENTITY trig_sufix "">` - change to the desired trigger suffix (relevant for 27 GeV)
- `<!ENTITY bes_phase "1">` - change to the desired BES phase (phase 2 prepared for but never run)
- `<!ENTITY dst "mu">` - change to the desired dst type (pico prepared for but never run)
- `<!ENTITY pions "false">` - change to true if you want to include pions in the tree
- `<!ENTITY out_path "/gpfs01/star/pwg/dneff/data/BES1/trees/">` - change to the desired output path **Most important**

In addition, all instances of `/star/u/dneff/git/tree_maker/` must be changed to the path of the user's repository.
All instances of `/star/u/dneff/gpfs/data/BES1/flow_flat_coefs/` must be changed to the path of the flow coefficients files.
Similarly, all instances of `/gpfs01/star/pwg/dneff/data/BES1/trees/` must be changed to the desired output path for logs/scripts/lists. This was typically the same as the output directory.
It was not possible to use xml variables in <Generator> tags, so these monotonous changes must be made manually.

Once the xml file is edited, submit the jobs to the condor queue with the following command:
```
star-submit <xml_file>
```

For 7 GeV Tree Maker, for example, from the subs/ directory:
```
star-submit sub7.xml
```

### Order of Running

Each step of this process must be run in order, as the output of the previous step is used as input for the next step.
The order of running is as follows:
1. Phi Flattening
2. Event Plane Flattening
3. Tree Making

### Output

The output root trees along with QA root files should be found in the output directory specified in the submission xml file.