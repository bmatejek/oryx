![alt text](oryx-logo.png)

Oryx is a Python package to extract the wiring diagram from segmented EM images with corresponding synapse labelings.
The function `ExtractWiringDiagram` takes as input two arguments: the segmentation and a synapse mask. The function returns the corresponding skeleton, distance from each synapse endpoint to the soma, and width of the neuron at each skeleton joint.