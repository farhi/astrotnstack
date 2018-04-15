ASTROTNSTACK: automatically align and stack astro-photography pictures
======================================================================

![Image of AstRotNStack](https://github.com/farhi/astrotnstack/blob/master/doc/screenshot.png)

**For a more recent implementation, refer to https://github.com/farhi/matlab-mastrostack**

This function (for Matlab/Octave) gets a list of images, and automatically determines bright stars as control points. These are followed along pictures, and used to build an affine transformation at constant scale (e.g. a rotation and translation). All images are stacked over the first image used as reference. The images can be given as file names (may include wildcards), or matrices from e.g. imread, and support both RGB and gray images.

As stars are used for the alignment, this method is suited for deep sky images, but not for planetary imaging.

This function does not make use of the phase correlation technique, but operates directly on the raw images. It assumes that at least two bright stars remain on each picture, but tolerate translations and rotations, as well as stars/control points that appear/disappear on the field of view. The procedure also ignores very sharp peaks, assuming they originate from dead pixels (and would then be static over images, leading to traces).

It is highly recommended to specify a 'dark' frame filename, which will be subtracted to all images, in order to e.g. remove most hot/dead pixels. To get such a 'dark', use your camera alone and shot once with the cap on, same settings as the other pictures (duration, ISO). Image should be full black.

You may as well specify a 'flat' frame filename, which will be divided to all images, in order to counter vignetting (darker borders). To get such a 'flat', shot once with the scope pointing at a uniform view (sky, white wall). Adapt the duration so that you get a rather gray image (not full black or white).

During the procedure, the current image and its stars/control points are indicated, as well as the stacked image. A waitbar indicates the progress. To abort the procedure, close the waitbar. The stacked image so far will be saved and returned. Supported image formats include JPG, PNG, TIFF, FITS.

The resulting stacked image is stored using the first/reference image name, adding '_stacked', in the PNG format. This routine can stack a hundred 4k images within 15 minutes. Memory requirements are about 2Gb for 4k images.

We recommend that you further use DarkTable or RawTherapee to enhance the low intensity features (black, brightness, contrast).

Syntax
------

```matlab
>> [stacked, output] = astrotnstack(images, dark, flat, N)
>> stacked = astrotnstack;
```

**input:**

- images: a filename name or cell of file names, that may include wildcards, or
      matrix (rgb, gray). When not given or empty, a file selector pops-up.
- dark: a single dark frame which is subtracted to all images.
    When not given or empty, a file selector pops-up. Use 'cancel' button
    to proceed without dark subtraction.
- flat: a single flat frame which is divided to all images. 
    When not given or empty, a file selector pops-up. Use 'cancel' button
    to proceed without flat normalisation.
- N:    the max number of control points to use. N=20 default.

Input can also be given as a structure with fields:

- N:         max number of control points              (default is 20)
- tol_rot:   rotation    tolerance in degrees       (default is 3 deg)
- tol_trans: translation tolerance     (default is 0.01 = 1% of width)
- test:      when 1, images are analysed, but result is not written.
- silent:    when 1, no image/wait bar is displayed (faster).
- dark:      indicates dark frame (single image)
- flat:      indicates flat frame (single image)
- images:    images to stack (filename, may use wildcard, or matrix)

or even as unsorted name/value pairs such as in:
```matlab
>> [stacked, output] = astrotnstack('images','*.JPG', 'dark','Dark.PNG', 'N',15)
```

**output:**

- stacked:  stacked images.
- output:   a structure with the file names, rotations and translations.

Examples:
---------

```matlab
>> astrotnstack('*.PNG', 'Dark.png', 'Flat.PNG');
>> astrotnstack({'*.PNG', imread('file.jpg') }, 'Dark.png', 'Flat.PNG');
>> astrotnstack('*.PNG', 'Dark.png', 'Flat.PNG', 'silent', 1);
```

Requirements/Installation
-------------------------

You only need Matlab, no external toolbox.

Copy the directory and navigate to it. Then type from the Matlab prompt:

  ```matlab
  addpath(pwd)
  ```
  
  If you also have **readraw** installed and available, you will be able to import
  RAW camera images.

Credits
-------
http://nghiaho.com/?page_id=671

See also: 

- LxNstack         https://sites.google.com/site/lxnstack/home
- Deep Sky Stacker http://deepskystacker.free.fr/french/
- Rot'n Stack      http://www.gdargaud.net/Hack/RotAndStack.html
- DarkTable        http://www.darktable.org/
- RawTherapee      http://rawtherapee.com/
- mastrostack      https://github.com/farhi/matlab-mastrostack
- readraw          https://github.com/farhi/matlab-readraw

E. Farhi Dec 2017. Version 1.6.2. GPL2.
