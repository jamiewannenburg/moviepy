Using the MultiCam class
------------------------

Multi Cam is when more than one camera records the same scene. For instance in a dialog scene you can switch between two cameras for the two people. It can also be used for a bullet time effect.

### Simple example

In this example two files from two cameras are called clip1.MP4 and clip2.MP4. (It is assumed that they started recoding at exactly the same moment.) The following script will start at 0s with clip1.MP4. At 5s it will switch to clip2.MP4. At 10s it will switch back to clip1.MP4 and end at 15s. The clips are concatenated and the returned clip is a VideoClip/CompositeClip as returned from the concatenate function.

.. literalinclude:: ../../examples/star_worms.py #include so files

    import MultiCam

    files = ['clip1.MP4','clip2.MP4']
    times = [
                [0,0],
                [5,1],
                [10,0],
                [15]
            ]
    clip = MultiCam.MultiCam(files,times)
    


### files

It often happens when recoding in a multi cam setup that your files are stored in the following way. There is a folder for each scene or day, under this folder there is a folder of each camera on set. Inside are files which may or may not be related to each other (shot at the same time) but are ordered in the naming convention of the camera.

This example shows how to use the get_files function in such a situation.

The function check_files, checks that the number of files in each folder is the same, then checks that the corresponding files have a similar file size. It then return a list of lists which can be used in the MultiCam function.

### sync

A common problem editors face is when they have files from multiple cameras which did not start at the same time. The audio from the cameras are often used to sync the files using plugins such as Plural Eyes. This function accomplishes the same thing. It calculates the shift in seconds relative the first file and stores it in the shift attribute.





