"""
This module is a collection of functions to help with multi cam editing, including the multicam class.

Note:
You need scipy for sync
"""
import os
import subprocess
#from moviepy.video.compositing.CompositeVideoClip import CompositeVideoClip
#from moviepy.video.VideoClip import VideoClip
from moviepy.audio.io.AudioFileClip import AudioFileClip
from moviepy.video.compositing.concatenate import concatenate
from moviepy.video.io.VideoFileClip import VideoFileClip
from moviepy.conf import FFMPEG_BINARY
import numpy as np
from scipy.io import wavfile
from scipy.signal import fftconvolve
from tempfile import mkdtemp
import shutil
import pygame as pg

class MultiCam:
    """
    This class
    
    :param clips: Either a list of VideoFileClips, or a list of filenames where clips are situated. 
    
    :param times: list of lists where the first value is the time and the second is the camera number. If left empty or None the get_times function should be called before get_clip.
    
    :param slowmo: dictionary where the keys are the index of the clip that is slowed down, the value is the amount it was slowed down.
    
    :param *kwargs: The rest of the parameters are the optional arguments for the concatinate function.
    """
    def __init__(self,clips,times=None,shift=[],slowmo={}):
        ### TODO add audio settings here????
        self.shift = shift
        self.times = times
        self.converted = False
        if isinstance(clips[0],VideoFileClip):
            self.clips = clips
            filenames = []
            for clip in clips:
                filenames.append(clip.parameters['filename'])
        else:
            # initialise clips later to save memory
            self.clips = None
            filenames = clips
        # seperate filename and extension
        self.filenames = []
        for f in filenames:
            f_a = f.split('.')
            if len(f_a) < 2:
                raise Exception("No extension?")
            if len(f_a) == 2:
                self.filenames.append([f_a[0],'.'+f_a[1]])
            else:
                f1 = []
                for fp in f_a[:-1]:
                    f1 = f1 + fp
                self.filenames.append([f1,'.'+f_a[-1]])
                
        self.slowmo = slowmo
        
    def load_clips(self):
        """ Load clips from filenames.
        """
        self.clips = []
        for f in self.filenames:
            self.clips.append( VideoFileClip(f[0]+f[1]) )
        
    
    def sync(self,fps=11025,nbytes=2,low_memory=False,print_progress=False, convert=False):
        """
        This function calculates the shift neccisary for the other cameras to be in sync
        with the first camera. It uses scipy's fftconvolve to compute the 
        cross correlation.
        :param convert: if convert is True, the audio from the video file is written to a wave file. (This uses scipy to read the file if it exists.)
        """
        # first file (refence)
        if convert:
            # only use wav if convert is on
            if os.path.exists(self.filenames[0][0]+'.wav'):
                with open(self.filenames[0][0]+'.wav','rb') as f:
                    fs,data = wavfile.read(f)
                # see if settings changed
                if fs != fps:
                    data = write_audio(self.filenames[0],fps,nbytes,overwrite=True)
            else:
                data = write_audio(self.filenames[0],fps,nbytes,overwrite=True)
                
        else:
            clip = AudioFileClip(self.filenames[0][0]+self.filenames[0][1])
            data = clip.to_soundarray(fps=fps, nbytes=nbytes)[0] #### is this right
            del clip.reader ############### maak seker
        
        if low_memory:
            reference = np.memmap(self.filenames[0][0]+'.dat', dtype='int16', mode='w+',shape=data.shape)
            reference = data[:]
            del data
        else:
            reference = data[:]
            del data
        
        # the rest (to sync)
        shift = [] 
        for i in range(len(self.filenames)-1):
            if print_progress:
                print "Syncing "+str(i+2)+" of "+str(len(self.filenames))
            
            
            if convert:
                # only use wav if convert is on
                if os.path.exists(self.filenames[i][0]+'.wav'):
                    with open(self.filenames[i][0]+'.wav','rb') as f:
                        fs,data = wavfile.read(f)
                    # see if settings changed
                    if fs != fps:
                        data = write_audio(self.filenames[i],fps,nbytes,overwrite=True)
                else:
                    data = write_audio(self.filenames[i],fps,nbytes,overwrite=True)
                    
            else:
                clip = AudioClip(self.filenames[i][0]+self.filenames[i][1])
                data = clip.to_soundarray(fps=fps, nbytes=nbytes)[0]
                del clip.reader
                
            if low_memory:
                to_sync = np.memmap(self.filenames[i][0]+'.dat', dtype='int16', mode='w+',shape=data.shape)
                to_sync = data[:]
                del data
            else:
                to_sync = data[:] ########### neccisary? (wrong)
                del data
            
            sync_time = get_shift(reference,to_sync,fps,low_memory=low_memory)
            
            if print_progress:
                print sync_time
            shift.append( sync_time )

        self.shift = shift
        return shift
    
    def get_clip(self,**kargs):
        """
        """
        if self.times == None:
            raise Exception('Times not specified. Run get_times.')
        if self.clips == None:
            self.load_clips()
        
        clip_array = []
        for i,time in enumerate(self.times[:-1]):
            
            clip_start = time[0]
            clip_end = self.times[i+1][0] 
            
            if i in self.slowmo:
                clip_start = clip_start*self.slowmo[i]
                clip_end = clip_end*self.slowmo[i]
            
            if time[1] > 0:
                clip_start = clip_start - self.shift[time[1]-1]
                clip_end = self.times[i+1][0] - self.shift[time[1]-1]
            # round frames?
            if clip_start < 0:
                clip_start = 0
            if clip_end < 0:
                clip_end = 0
            if clip_start > self.clips[time[1]].duration:
                clip_start = self.clips[time[1]].duration - 1.0/self.clips[time[1]].fps
            if clip_end > self.clips[time[1]].duration:
                clip_end = self.clips[time[1]].duration
                
            clip_array.append(self.clips[time[1]].subclip(clip_start,clip_end))
            
        return concatenate(clip_array,**kargs)
        
    def get_times(self):
        """
        This is not perfect but wil help you. Consider running this from the terminal.
        """
        self.load_clips()
        clip = self.clips[0] # reference ???????????????
            
        print "Click when you want to switch views."
        times = clip.preview(fps=5,audio=False,func=get_time)
        print "Now pick viewpoints."
        #TODO print "Use up and down arrow to select view points, enter to select"
        for i in len(times):
            clip.show(times[i][0])
            in_str = raw_input( "To which camera should it switch here?" )
            times[i].append(int(in_str))
            
        self.times = times
        return times
        
# is hierdie nodig? Hy laai hom inelkgeval in memory...
def write_audio(filename,fps,nbytes,overwrite=False):
    """
    :param filename: list with first index the filename, second the extension.
    """
    input_str = 'n\n'
    if overwrite:
        input_str = 'y\n'
            
    cmd = [ FFMPEG_BINARY, '-i', filename[0]+filename[1], '-vn',
                       '-acodec', 'pcm_s%dle'%(8*nbytes),
                       '-ar', '%d'%fps,
                       '-ac', '1', filename[0]+'.wav']
                
    p = subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    ret = p.communicate(input=input_str)
    
    with open(filename[0]+'.wav','rb') as f:
        fs,data = wavfile.read(f)
        
    return data
    
def check_numbers(base_folder,cameras,extension = '.MP4',raise_error=True):
    """
    Checks number of files in each camera folder. 
    If they differ raise an error with helpfull 
    information if raise_error is True.
    Otherwise simply prints info and returns False.
    """
    structure = get_file_structure(base_folder,cameras,extension=extension)
    
    num_of_files = len(structure[0])
    for i,gp in enumerate(structure):
        if len(gp) != num_of_files:
            if raise_error:
                raise Exception("The " + str(i+1) + "th folder has " + str(len(gp)) + " video files, while the previous has " + str(len(structure[i-1])))
            else:
                print "The " + str(i+1) + "th folder has " + str(len(gp)) + " video files, while the previous has " + str(len(structure[i-1]))
                return False
            
    return True
    
def check_sizes(base_folder,cameras,extension = '.MP4', tollerance=0.7, raise_error=True):
    """
    This function check that corresponding file sizes in camera folders are
    within a certain tollerance.
    If rease_error is True if fails with an error,
    if it is False it goes into an interactive mode and returns False only if a
    difference in filesize is not acceptible.
    """
    filenames = get_file_structure(base_folder,cameras,extension=extension)
    
    for i in range(len(filenames[0])):
        filesize_base = os.path.getsize(filenames[0][i])
        for j in range(len(filenames)-1):
            filesize = os.path.getsize(filenames[j+1][i])
            ok = False
            if float(min(filesize_base,filesize))/max(filesize_base,filesize) > tollerance:
                ok = True
            elif float(min(filesize_base,2*filesize))/max(filesize_base,2*filesize) > tollerance:
                # allow double frame rate on camera 1
                ok = True
            elif float(min(2*filesize_base,filesize))/max(2*filesize_base,filesize) > tollerance:
                # allow double frame rate on other cameras
                ok = True
            elif float(min(filesize_base,4*filesize))/max(filesize_base,4*filesize) > tollerance:
                # allow 4x frame rate on camera 1
                ok = True
            elif float(min(4*filesize_base,filesize))/max(4*filesize_base,filesize) > tollerance:
                # allow 4x frame rate on other cameras
                ok = True
            
            if not ok:
                print "File sizes to not match:"
                if raise_errors:
                    raise Exception(filenames[0][i]+ self.videoextension + " is " + str(filesize_base) + " but " + filenames[j+1][i]+ self.videoextension + " is " + str(filesize))
                else:
                    print filenames[0][i]+ self.videoextension + " is " + str(filesize_base) + " but " + filenames[j+1][i]+ self.videoextension + " is " + str(filesize)
                    cont = True
                    while cont:
                        ret = raw_input("Is this ok? [Y/n] ")
                        if ret == '' or ret == 'y' or ret == 'Y':
                            cont = False
                        elif ret == 'n' or ret == 'N':
                            return False
    return True
    
#    def __del__(self):
#        if self.clips != None:
#            for clip in self.clips:
#                clip.reader.close()

def check_files(base_folder,cameras,extension = '.MP4',raise_error=True):
    """
    Used to make sure video files in corresponding 
    camera folders correspond to one another.
    Runs check_numbers and check_sizes.
    """
    if not check_numbers(base_folder,cameras,extension=extension,raise_error=raise_error):
        return False
        
    if not check_sizes(base_folder,cameras,extension=extension,raise_error=raise_error):
        return False
        
    return True

def get_file_structure(base_folder,cameras,extension='.MP4'):
    """
    Returns directory listing matching extension for each camera.
    """
    structure = []
    for i,camera in enumerate(cameras):
        path = os.path.abspath(os.path.join(base_folder, camera)) 
        structure.append( [] )
        for f in sorted(os.listdir(path)):
            if extension in f:
                clip_path = os.path.join(path,f)
                structure[i].append(clip_path)
    return structure
    
def get_files(base_folder,cameras,extension='.MP4'):
    """
    Writes and returns the directory structure (without file extension).
    The file structure is a list whose indecies correspond to the cameras' indecies.
    Each index of which can be used to initiate a MultiCam class.
    eg. 
    [
        ['/basefolder/camera1/GOPR01234.MP4','/basefolder/camera2/GOPR01239.MP4'],
        ['/basefolder/camera1/GOPR01235.MP4','/basefolder/camera2/GOPR01241.MP4']
    ]
    Note: to make sure files are consistent run check_files before this function.
    """
    # first get the directory listing of each camera
    structure = get_file_structure(base_folder,cameras,extension=extension)
                
    # now invert to get correct structure
    filenames = []
    for i in range(len(structure[0])):
        filenames.append([])
        for j in range(len(cameras)):
            filenames[i].append(structure[j][i])
        
    return filenames
    

def get_shift(reference,to_sync,sr,low_memory=False,max_shift=0.33333):
    """
    returns the shift neccisary for to_sync to be in sync with reference. 
    :param reference: 
    :param to_sync: 
    :param sr: 
    :param low_memory: If low_memory is True, the computations will be done with numpy.mmap arrays.
    :param max_shift: 
    """
    s1 = reference.shape[0]
    s2 = to_sync.shape[0]
    mylength = (s1+s2-1)
    for i in range(1000):
        if mylength<2**i:
            mylength = 2**(i-1)
            break
    try:
        if low_memory:
            tmp_dir = mkdtemp()
            l1f = path.join(tmp_dir, 'l1.dat')
            l2f = path.join(tmp_dir, 'l2.dat')
            l1 = np.memmap(l1f, dtype='int16', mode='w+',shape=(mylength))
            np.put(l1,range(mylength),reference)
    
            l2 = np.memmap(l2f, dtype='int16', mode='w+',shape=(mylength))
            np.put(l2,range(mylength),to_sync[:])
        else:
            l1 = np.zeros(mylength)
            l1[:len(reference)] = reference
            l2 = np.zeros(mylength)
            l2[:len(to_sync)] = to_sync
            
        r = fftconvolve(l1,l2[::-1])
        maxplace_rel = np.argmax(abs(r[mylength-int(sr*max_shift):mylength+int(sr*max_shift)]))
        maxplace = maxplace_rel + mylength - int(sr*max_shift)
        
        sync_time = float(maxplace-mylength)/float(sr) 
        
    finally:
        if low_memory:
            try:
                shutil.rmtree(tmp_dir)  # delete directory
            except OSError as exc:
                if exc.errno != 2:  # code 2 - no such file or directory
                    raise  # re-raise exception

    return sync_time

def get_time(event,img,t):
    """
    Returns time clicked
    """
    if event.type == pg.MOUSEBUTTONDOWN:
        x,y = pg.mouse.get_pos()
        print t
        return [t]

