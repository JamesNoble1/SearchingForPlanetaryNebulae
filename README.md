# SearchingForPlanetaryNebulae
Code written for University Project entitled 'Searching for Planetary Nebulae and testing the application of machine learning'

REMOVE_WEIGHTMAP.py changes the download script from ESO website so only necessary files are downloaded.

ReOrg4.py Matches the downloaded files in to H alpha and red pairs (showing the same area of sky) and send them to a subtraction program.

VST_cutout_all.py creates images of PN from the subtracted images.  These images are used to train the ML model.

Object_detection.py loads a saved model and checks input images for PN.
