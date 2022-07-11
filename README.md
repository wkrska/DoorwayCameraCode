# DoorwayCameraCode
This repo contains everything you need to replicate the results found in the main document and supplemental document. 
[doorway_camera.m](./doorway_camera.m) is the main function that will perform the reconstrution. It requires several parameters to be passed into it in order to perform the reconstrution. 

[run_doorway_camera.m](./run_doorway_camera.m) has all the function calls needed to run each reconstruction shown in the main document and supplemenary document. This file also explains in detail what each parameter is. NOTE: Each reconstruction takes about 15 minutes to run on a moderately powerful consumer CPU given our chosen parameters. 

[run_plot_pdfs.m](run_plot_pdfs.m) will create the reconstruction images seen in both documents. The specific parameters of the reconstruction you are plotting must be specified in the file, and this reconstruction must have already been performed, as each reconstruction saves a unique datafile. 

To run reconstructions on new data, simply create a new folder in the same format as one of the examples, and add the filepath to the switch-block(s) in [doorway_camera.m](./doorway_camera.m)
