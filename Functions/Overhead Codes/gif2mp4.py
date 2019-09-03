import moviepy.editor as mp
from os import listdir
from os import remove

def convert(fname):
    clip = mp.VideoFileClip(fname)
    clip.write_videofile(fname[:-3]+'mp4')
    remove(fname)

if __name__ == "__main__":
    afiles = listdir()
    for i in afiles:
        if i.find('.gif') is not -1:
            convert(i)

