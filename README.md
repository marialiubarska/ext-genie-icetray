# external-genie-icetray

**structure:**

```
genie-generator
      GenieDriversIceCube
      HelperClasses 
      genie
      ...
hepevt-reader
      ...
```

**TODO:**

/srv -- why does it look for LHAPDF there??

add cross-sections key to gicegen ($GSPLOAD was removed between R-2_6_0 and R-2_12_8)

add exeptions class for reading (like in R-2_6_0) - do I need this? (anyway, redo args reading to look nicer) 

finish changes for R-2_12_8 (+ plus try to install my R-2_12_8?)

install R-3_X_X (latest or corresponding to GH? are they the same?)

change ext-genie to be able to use it with R-3_X_X (can you switch models in it?)

install GENIE-HEDIS (or use existing installation?? -- then need to be careful with externals)

adapt changes for R-3_X_X to work with GH

...


-----

will it actually work without container??

*how to make code GENIE version independent?*
