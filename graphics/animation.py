import imageio



def Images2Gif(filename, images, duration=0.5):
    imageio.mimsave(filename, images, duration=duration)
