function [im0, output] = astrotnstack(varargin)
% ASTROTNSTACK: automatically align and stack astro-photography pictures.
%
% This function gets a list of images, and automatically determines bright stars
%  as control points. These are followed along pictures, and used to build an 
%  affine transformation at constant scale (e.g. a rotation and translation).
%  All images are stacked over the first image used as reference. The images can
%  be given as file names (may include wildcards), or matrices from e.g. imread,
%  and support both RGB and gray images.
%
% As stars are used for the alignment, this method is suited for deep sky 
%  images, but not for planetary imaging.
%
% This function does not make use of the phase correlation technique, but
%  operates directly on the raw images. It assumes that at least two bright  
%  stars remain on each picture, but tolerate translations and rotations, as
%  well as stars/control points that appear/disappear on the field of view.
% The procedure also ignores very sharp peaks, assuming they originate from dead
%  pixels (and would then be static over images, leading to traces).
%
% It is highly recommended to specify a 'dark' frame filename, which will be
%  subtracted to all images, in order to e.g. remove most hot/dead pixels.
%  To get such a 'dark', use your camera alone and shot once with the cap on, same
%  settings as the other pictures (duration, ISO). Image should be full black.
%
% You may as well specify a 'flat' frame filename, which will be divided
%  to all images, in order to counter vignetting (darker borders).
%  To get such a 'flat', shot once with the scope pointing at a uniform view 
%  (sky, white wall). Adapt the duration so that you get a rather gray image 
%  (not full black or white).
%
% During the procedure, the current image and its stars/control points are 
%  indicated, as well as the stacked image. A waitbar indicates the progress.
% To abort the procedure, close the waitbar. The stacked image so far will be
%  saved and returned.
% Supported image formats include JPG, PNG, TIFF, FITS.
%
% The resulting stacked image is stored using the first/reference image name,
%  adding '_stacked', in the PNG format. This routine can stack a hundred 4k
%  images within 15 minutes. Memory requirements are about 2Gb for 4k images.
%
% We recommend that you further use DarkTable or RawTherapee to enhance the low
%  intensity features (black, brightness, contrast).
%
% Syntax: [stacked, output] = astrotnstack(images, dark, flat, N)
%         stacked = astrotnstack;
%
% input:
%   images: a filename name or cell of file names, that may include wildcards, or
%           matrix (rgb, gray). When not given or empty, a file selector pops-up.
%   dark: a single dark frame which is subtracted to all images.
%         When not given or empty, a file selector pops-up. Use 'cancel' button
%         to proceed without dark subtraction.
%   flat: a single flat frame which is divided to all images. 
%         When not given or empty, a file selector pops-up. Use 'cancel' button
%         to proceed without flat normalisation.
%   N:    the max number of control points to use. N=20 default.
%
%   Input can also be given as a structure with fields:
%           N:         max number of control points              (default is 20)
%           tol_rot:   rotation    tolerance in degrees       (default is 3 deg)
%           tol_trans: translation tolerance     (default is 0.01 = 1% of width)
%           test:      when 1, images are analysed, but result is not written.
%           silent:    when 1, no image/wait bar is displayed (faster).
%           dark:      indicates dark frame (single image)
%           flat:      indicates flat frame (single image)
%           images:    images to stack (filename, may use wildcard, or matrix)
%
%   or even as unsorted name/value pairs such as in:
%    [stacked, output] = astrotnstack('images','*.JPG', 'dark','Dark.PNG', 'N',15)
%
% output:
%   stacked:  stacked images.
%   output:   a structure with the file names, rotations and translations.
%
% Example:
%   astrotnstack('*.PNG', 'Dark.png', 'Flat.PNG');
%   astrotnstack({'*.PNG', imread('file.jpg') }, 'Dark.png', 'Flat.PNG');
%   astrotnstack('*.PNG', 'Dark.png', 'Flat.PNG', 'silent', 1);
%
% Credits: http://nghiaho.com/?page_id=671
%
% See also: 
%           LxNstack         https://sites.google.com/site/lxnstack/home
%           Deep Sky Stacker http://deepskystacker.free.fr/french/
%           Rot'n Stack      http://www.gdargaud.net/Hack/RotAndStack.html
%           DarkTable        http://www.darktable.org/
%           RawTherapee      http://rawtherapee.com/
%
% E. Farhi Dec 2017 GPL2. Version 1.6.2

% init options
options.tol_trans = 0.01; % 1% of total image width
options.tol_rot   = 3;    % in degrees
options.N         = 20;
options.dark      = '';
options.flat      = '';
options.test      = 0;
options.images    = '';
options.silent    = 0;
fields = fieldnames(options);
f1 = '';
f2 = '';

% get input arguments ----------------------------------------------------------
index = 1;
while index <= numel(varargin)
  % check for numeric -> N
  % check for option as name/value pair
  % check as struct -> options
  % otherwise store in order dark, flat
  this = varargin{index};
  if isnumeric(this) && isscalar(this)
    N = this;
  elseif ischar(this) && any(strcmp(this, fields)) && index < numel(varargin)
    val = varargin{index+1};
    f   = fields{strcmp(this, fields)};
    options.(f) = val;
    index=index+1;
  elseif (ischar(this) || iscellstr(this))
    if     isempty(options.images),   options.images    = this;
    elseif isempty(options.dark) options.dark = this;
    elseif isempty(options.flat) options.flat = this;
    end
  elseif isstruct(this)
    for f=fieldnames(this)'
      options.(f{1}) = this.(f{1});
    end
  end
  index = index+1;
end

% pass 'options' into names variables
s        =options.images;
N        =options.N; 
tol_trans=options.tol_trans; 
tol_rot  =options.tol_rot;
dark     =options.dark; 
flat     =options.flat;
test     =options.test; 
im0      = []; output = [];

% get file names ---------------------------------------------------------------
disp([ mfilename ': E. Farhi (2017) GPL2']);

dark = resolvefiles(dark, 'Pick a Dark Frame or Cancel to Ignore');
if ~isempty(dark) && iscellstr(dark) && ~isequal(dark, 0)
  disp([ mfilename ': subtracting dark frame ' dark{1} ' for all images.' ]);
  dark = imread(dark{1});
else
  disp([ mfilename ': WARNING: not using any dark frame.' ])
  disp('    TIP: Take a shot with cap on and use it to remove e.g. dead pixels.');
  dark = [];
end

flat = resolvefiles(flat, 'Pick a Flat Frame or Cancel to Ignore');
if ~isempty(flat) && iscellstr(flat) && ~isequal(flat, 0)
  disp([ mfilename ': dividing flat frame ' flat{1} ' for all images.' ]);
  flat = imread(flat{1});
  flat = double(rgb2gray(flat));
  flat = flat/max(flat(:)); % normalise to max
else
  disp([ mfilename ': WARNING: not using any flat frame.' ])
  disp('    TIP: Take a shot on a white surface or uniform sky to compensate');
  disp('         for vignetting (darker on borders).');
  flat = [];
end

s      = resolvefiles(s, 'Pick images to stack');  % file names
if isempty(s) || isequal(s, 0) || ~iscell(s)
  disp('Nothing to do. Exiting');
  return;
end

% initiate treatment -----------------------------------------------------------
disp([ mfilename ': will treat ' num2str(numel(s)) ' files.' ]);

points.x = [];  % location  of 'stars', size [ numel(s) N ]
points.y = [];  % location  of 'stars', size [ numel(s) N ]
points.m = [];  % intensity of 'stars', size [ numel(s) N ]
M0       = 1;

if ~options.silent
  f1  = figure('Name','Current image');
  f2  = figure('Name','Stacked image');
  delete(findall(0, 'Tag', [ mfilename '_waitbar' ]));
  wb  = waitbar(0, 'Stacking images...'); set(wb, 'Tag', [ mfilename '_waitbar' ]);
end
im0 = 0;
im  = [];
output.angle  = [ 0 ];
output.T      = [ 0 0 ];
output.files  = {};
output.R      = {};
output.points = {};
output.time   = [];

t0 = clock;
totaltexposure = 0;
totalimages    = 0;
reference      = '';

for index=1:numel(s)

  if ~options.silent
    if ~ishandle(wb)
      disp('Aborting.')
      break;
    end
    try
      waitbar(index/numel(s), wb, s{index});
    end
  end
  
  % ========== read image
  if ischar(s{index})
    name = s{index};
    try
      im = imread(name);
    catch
      disp([ mfilename ': Ignoring ' name ' (not a supported image format).' ])
      continue;
    end
  elseif isnumeric(s{index})
    im   = s{index};
    name = sprintf('data_%i', totalimages + 1);
  else
    continue
  end
  if ~isinteger(im) && isreal(im) || ~isa(im, 'uint8')
    im = uint8(double(im)*256/max(double(im(:))));  % make it a uint8
  end
  totalimages = totalimages + 1;
  
  % remove dark hot spots
  if ~isempty(dark)
    im=im-dark; 
  end
  
  % compensate vignetting (borders)
  if ~isempty(flat)
    for l=1:size(im,3)
      layer = double(im(:,:,l));
      layer = layer./flat;
      im(:,:,l) = cast(layer, class(im));
    end
  end
  
  if totalimages == 1 % reference
    disp([ mfilename ': ' num2str(index) '/' num2str(numel(s)) ': ' name ]);
    meta = display_metainfo(s{index});
    if isempty(meta), continue; end
    im0  = uint16(im);
    if isfield(meta, 'ExposureTime')
      totaltexposure = totaltexposure+meta.ExposureTime; 
      output.time = [ output.time meta.ExposureTime ];
    end
    output.files{end+1} = s{index};
    reference           = name;
  else
    if totalimages > 2
      dt_from0     = etime(clock, t0);
      dt_per_image = dt_from0/(index-2);
      % remaining images: numel(s)-index
      eta    = dt_per_image*(numel(s)-index+1);
      ending = addtodate(now, ceil(eta), 'second');
      ending = [ 'Ending ' datestr(ending) ];
      eta    = sprintf('ETA %i [s]. %s', round(eta), ending);
    else eta=''; ending='';
    end
    disp([ mfilename ': ' num2str(index) '/' num2str(numel(s)) ': ' name '. ' eta]);
    if ~options.silent
      if ishandle(wb)
        set(wb, 'Name', [ num2str(index) '/' num2str(numel(s)) ' ' ending ])
      else break;
      end
    end
    meta = display_metainfo(s{index}, false);
    if isempty(meta), continue; end
  end
  
  % ========== find control points on gray image
  points = find_control_points(im, totalimages, N, points, tol_trans*size(im,1));

  if totalimages == 1
    continue; 
  end
  
  % now we identify the points which keep distances and angles continuous
  % from one image to the next

  [x1,y1,x2,y2, p1_orig,p2_orig, p1_axis,p2_axis] = analyse_dist_angles(...
    points, totalimages, tol_trans*size(im,1), tol_rot);

  % get the best similarity
  if ~p1_orig || ~p2_orig || ~p1_axis || ~p2_axis
    disp([ mfilename ': WARNING: not enough control points for ' num2str(index) ' (axis). Will ignore image ' name ])
    continue;
  end
  
  % identify the current points which are common to the reference, and match
  % distances AND rotations
  p1 = p1_orig; p2=p2_orig;
  d1  = sqrt( (x1-x1(p1)).^2 + (y1-y1(p1)).^2 );
  t1  = atan2( y1-y1(p1),       x1-x1(p1))*180/pi;
  d2  = sqrt( (x2-x2(p2)).^2 + (y2-y2(p2)).^2 );
  t2  = atan2( y2-y2(p2),       x2-x2(p2))*180/pi;
  [ok1,ok2] = find_similar2(d1,             d2,             tol_trans*size(im,1), ...
                            t1-t1(p1_axis), t2-t2(p2_axis), tol_rot*pi/180);
  
  % plot image
  if ~options.silent && ishandle(f1), 
    figure(f1);
    [h,colors]=plot_im_points(im, name, ...
      points.x(totalimages,:), points.y(totalimages,:), points.m(totalimages,:));
    if all(ishandle(h))
      % we highlight p2_orig and p2_axis control points
      set(h(p2_orig),'MarkerFaceColor',colors(p2_orig));
      set(h(p2_axis),'MarkerFaceColor',colors(p2_axis));
      drawnow
    end
  end
  
  if numel(ok1) <= 1 || numel(ok2) <= 1
    disp([ mfilename ': WARNING: not enough control points for ' num2str(index) ' (common). Will ignore image ' name ])
    continue
  end
  
  % we make a check for wrong guesses
  delta = (t1(ok1)-t1(p1_axis)) - (t2(ok2)-t2(p2_axis));
  bad   = find(abs(delta) > tol_rot);
  ok1(bad) = [];
  ok2(bad) = [];
  if numel(ok1) <= 1 || numel(ok2) <= 1
    disp([ mfilename ': WARNING: not enough control points for ' num2str(index) ' (common2). Will ignore image ' name ])
    continue
  end
  
  % compute the mean translation (x,y) and rotation angle wrt to reference image
  x1  = x1(ok1); y1=y1(ok1);
  x2  = x2(ok2); y2=y2(ok2);
  
  % theta2-theta1 is an estimate of the rotation angle
  % compute the affine transformatin
  [ret_R, ret_t] = rigid_transform_3D([x2 ; y2]', [x1 ; y1]');
  if isempty(ret_R)
    disp([ mfilename ': WARNING: invalid affine transformation. Skipping.'])
    continue
  end
  % compute an estimate of the translation and rotation from the identified
  % control points orig and axis
  theta1 = t1(p1_axis); % t1(p1_orig) == 0
  theta2 = t2(p2_axis);
  theta = asind(ret_R(1,2));
  if abs(theta-(theta2-theta1)) > tol_rot/3
    disp([ mfilename ': WARNING: invalid affine rotation. Skipping.']);
    continue
  end
  
  output.angle(end+1) = theta;
  output.T(end+1,1:2) = ret_t;
  output.R{end+1}     = ret_R;
  output.files{end+1} = name;
  output.points{end+1}= points;
  disp([ '  Angle=' num2str(theta,3) ' [deg] ; Translation=' mat2str(ret_t,3) ' [pix]']);
  
  % apply the rotation and stack raw counts
  if isempty(test) || test==0
    [im,M] = imaffine(im, ret_R, ret_t);
    M0=M0+M;
    clear M
    im0 = im0 + uint16(im); % add rotated/translated image on the reference
  end
  if isfield(meta, 'ExposureTime')
    totaltexposure = totaltexposure+meta.ExposureTime; 
    output.time    = [ output.time meta.ExposureTime ];
  end
  clear im
  
  % evaluate the stacked image and plot it
  if ~options.silent && ishandle(f2)
    figure(f2);
    if ischar(reference)
      [p,f] = fileparts(reference);
    else f=''; end
    [h,colors]=plot_im_points(eval_stacked(im0, M0), [ 'Stacked ' f ], ...
      points.x(1,:), points.y(1,:), points.m(1,:));
    if all(ishandle(h))
      set(h(p1_orig),'MarkerFaceColor',colors(p1_orig));
      set(h(p1_axis),'MarkerFaceColor',colors(p1_axis));
      drawnow
    end
  end
  
end % for index

if  ~options.silent && ishandle(f1), close(f1); end
if  ~options.silent && ishandle(wb), close(wb); end

meta.ExposureTime   = totaltexposure;
output.ExposureTime = totaltexposure;
disp([ mfilename ': Elapsed time ' num2str(etime(clock, t0)) ' [s]' ])
disp([ mfilename ': Total exposure on stacked image: ' num2str(totaltexposure) ' [s]' ]);

if isempty(test) || test==0
  % normalise to counts, enhance contrast
  im0 = eval_stacked(im0*ceil(2^16/max(im0(:))), M0);
  im0 = im0*ceil(2^16/max(im0(:)));
  f   = write_stacked( im0, reference, meta);
  if ~options.silent && ishandle(f2) 
    figure(f2); clf;
    set(f2, 'Name', [ 'log scale ' f ]);
    image(imlogscale(im0)); 
    title([ f ]);
    axes('position',[ 0.7 0.05 0.2 0.2 ]);
    h=imhist(im0); bar(h); title('Brightness histogram');
    set(gca, 'YScale','log');
  end
end

if ~options.silent
  f3=figure('Name',[ mfilename ': metrics' ]);
  plot_metrics(output);
end

% ------------------------------------------------------------------------------
% private function 
% ------------------------------------------------------------------------------

function [x1,y1,x2,y2, p1_orig,p2_orig, p1_axis,p2_axis] = analyse_dist_angles(points, index, tol_trans, tol_rot)

  % there is a global translation, and a minor rotation which brings some of
  % the current image back onto some of the previous one.
  x2 = points.x(index,:); x1=points.x(1,:);
  y2 = points.y(index,:); y1=points.y(1,:);
  i1 = find(x1>0 & y1>0);
  i2 = find(x2>0 & y2>0);
  x1 = x1(i1);  x2 = x2(i2); 
  y1 = y1(i1);  y2 = y2(i2);
  
  % we search for similar points, taking as reference one of the control points
  % the distances are p1_orig || ~p2_orig || ~p1_axis || ~p2_axis preserved in all images, and independent of rotations.
  ok1_best = 1; ok2_best = 1; p1_orig = 0; p2_orig = 0; p1_axis  = 0;  p2_axis = 0;
  for p1=1:numel(x1)
    % we select an initial control point in the reference image, that is used
    % as origin for the distances.
    % the reference can change when the second image moves and the overlap can
    % exclude some starts from the intersection.
    d1  = sqrt( (x1-x1(p1)).^2 + (y1-y1(p1)).^2 );
    
    for p2=1:numel(x2)
      % then we do the same on the current image
      d2  = sqrt( (x2-x2(p2)).^2 + (y2-y2(p2)).^2 );
      % d1(ok1) = d2(ok2)
      [ok1,ok2] = find_similar(d1, d2, tol_trans);
      
      if numel(ok1) > ok1_best && numel(ok2) > ok2_best 
        % store this solution as best bet, with at least 2 similar points
        ok1_best = numel(ok1); ok2_best = numel(ok2);
        p1_orig  = p1;         p2_orig  = p2;
      end
    end
  end
  
  if ~p1_orig || ~p2_orig, return; end
  
  % now we search for a second control point from which angles are measured
  
  ok1_best = 1; ok2_best = 1; 
  % compute absolute angle values using the same 'origin' as that for distances
  t1  = atan2( y1-y1(p1_orig), x1-x1(p1_orig))*180/pi;  
  t2  = atan2( y2-y2(p2_orig), x2-x2(p2_orig))*180/pi;
  for p1=1:numel(x1)
    if p1 == p1_orig, continue; end
    for p2=1:numel(x2)
      if p2 == p2_orig, continue; end
      % we try combinations so that angle(axis-orig) is used as reference
      % and other stars angles are measured from that direction.
      % this is then independent of the image rotations.
      [ok1,ok2] = find_similar(t1-t1(p1), t2-t2(p2), tol_rot*pi/180);
      % the match for distances and angle should coincide
      if numel(ok1) > ok1_best && numel(ok2) > ok2_best 
        ok1_best = numel(ok1); ok2_best = numel(ok2);
        p1_axis  = p1;         p2_axis  = p2;
      end
    end
  end
  
function info = display_metainfo(this, flag)
    if nargin <2, flag=true; end
    info=[]; 
    if ~ischar(this), return; end
    % print information
    try
      info  = imfinfo(this);
    catch
      return
    end
    clear this
    f = {'Software', 'DateTime', 'WhiteBalance', 'ExposureProgram', ...
      'LightSource', 'ExposureMode', 'ExposureBiasValue', 'FNumber', ...
      'ExposureTime', 'DigitalCamera', 'XResolution', 'YResolution', ...
      'ISOSpeedRatings','Author', 'Comment', ...
      'Make','Model','Source'}; 
    for i=1:numel(f)
      if isfield(info, 'DigitalCamera') && isfield(info.DigitalCamera, f{i})
        info.(f{i}) = info.DigitalCamera.(f{i});
      end
      if flag && isfield(info, f{i})
        try; disp([ '  ' f{i} ' = ' num2str(info.(f{i})) ]); end
      end
        
      if ~isfield(info, f{i}), info.(f{i})=''; end
    end
    
function im2 = eval_stacked(im2, M0)
  if isempty(im2), return; end
  M0  = uint16(M0);
  for z=1:size(im2,3)
    im2(:,:,z) = im2(:,:,z)./M0;
  end
  
function points = find_control_points(im, index, N, points, tol_trans)
  % find control points on gray image
  
  im= rgb2gray(im);
  points.x(index,:) = zeros(1, N);
  points.y(index,:) = zeros(1, N);
  points.m(index,:) = zeros(1, N);
  points.d(index,:) = zeros(1, N);
  points.t(index,:) = zeros(1, N);
  points.sx(index,:)= zeros(1, N);
  points.sy(index,:)= zeros(1, N);
  
  % find 'max' intensity locations. Every time we have a spot, we 'blank' it
  % around so that other ones are separated by 'tol_trans*10'
  for p=1:N
    [x1, y1, m1, im, sx, sy] = max_and_zero(im, tol_trans*10, tol_trans*10);
    if ~x1 || ~y1, continue; end % can not find a new star
    points.x(index,p) = x1;
    points.y(index,p) = y1;
    points.m(index,p) = m1;
    points.sx(index,p)= sx;
    points.sy(index,p)= sy;
  end

  % sort points with increasing 'y' value
  [points.y(index,:),p] = sort(points.y(index,:));
  points.x(index,:) = points.x(index,p);
  points.m(index,:) = points.m(index,p);
  points.sx(index,:)= points.sx(index,p);
  points.sy(index,:)= points.sy(index,p);
  
function [index_A_in_B, index_B_in_A] = find_similar2(A, B, find_tol1, C, D, find_tol2)
  % find_similar2: find indices in 'A' and 'C' that are similar to 'B' and 'D'
  %   within 'find_tol'
  %
  %   A(index_A_in_B) ~ B
  %   B(index_B_in_A) ~ A
  %
  % input:
  %   A:        vector for which values must be compared with the reference
  %   B:        list of 'reference' values
  %   find_tol1: tolerance used for comparison between A and B
  %   C:        vector for which values must be compared with the reference
  %   D:        list of 'reference' values
  %   find_tol2: tolerance used for comparison between C and D
  %
  % output:
  %   index_A_in_B: 'A' elements which are in 'B' within 'find_tol'
  %   index_B_in_A: 'B' elements which are in 'A' within 'find_tol'

  iA = zeros(size(A)); iB = zeros(size(B));
  index_A_in_B = [];
  index_B_in_A = [];
  
  for indexB=1:numel(B)
    x = B(indexB);
    [mindiff1, fA] = min(abs(A-x));
    y = D(indexB);
    [mindiff2, fC] = min(abs(C-y));
    if mindiff1 < find_tol1 && mindiff2 < find_tol2 && ~iA(fA) && ~iB(indexB)
      iA(fA) = 1; iB(indexB)=1;
      index_A_in_B = [ index_A_in_B fA ];
      index_B_in_A = [ index_B_in_A indexB ];
    end
  end
  
function [index_A_in_B, index_B_in_A] = find_similar(A, B, find_tol)
  % find_similar: find indices in 'A' that are similar to 'B' 
  %   within 'find_tol'
  %
  %   A(index_A_in_B) ~ B
  %   B(index_B_in_A) ~ A
  %
  % input:
  %   A:        vector for which values must be compared with the reference
  %   B:        list of 'reference' values
  %   find_tol: tolerance used for comparison
  %
  % output:
  %   index_A_in_B: 'A' elements which are in 'B' within 'find_tol'
  %   index_B_in_A: 'B' elements which are in 'A' within 'find_tol'
  
  iA = zeros(size(A)); iB = zeros(size(B));
  index_A_in_B = [];
  index_B_in_A = [];
  
  for indexB=1:numel(B)
    x = B(indexB);
    [mindiff, fA] = min(abs(A-x));
    if mindiff < find_tol && ~iA(fA) && ~iB(indexB)
      iA(fA) = 1; iB(indexB)=1;
      index_A_in_B = [ index_A_in_B fA ];
      index_B_in_A = [ index_B_in_A indexB ];
    end
  end
  
function [Y,M] = imaffine(X, ret_R, ret_t)
  % apply an affine transform to an image
  %
  % input:
  %  im: initial image
  %  R,t: rotation matrix(2x2) and translation
  %
  % output:
  %  Y: transformed image
  %  M: mask indicating assigned pixels in new image
	m = size(X,1);
	n = size(X,2);
	s = size(X); if numel(s) < 3, s(3)=1; end
	x1 = 1:m;
	y1 = 1:n;
	[x1,y1]= meshgrid(x1,y1);
	A = [ x1(:) y1(:) ];
	B = (ret_R*A') + repmat(ret_t, 1, size(A,1));
	clear A
	B = round(B); x2=B(1,:); y2=B(2,:); 
	clear B
	Y = X*0;
	M=uint8(zeros(m,n));
	x2(x2<1)= 1; y2(y2<1)=1;
	x2(x2>size(X,1)) =1;
	y2(y2>size(X,2)) =1;
	for z=1:s(3)
	  i1 = sub2ind(size(X),x1,y1,z*ones(size(x1)));
	  i2 = sub2ind(size(Y),x2,y2,z*ones(size(x2)));
	  Y(i2) = X(i1);
	  if z == 1, M(i2) = 1; end
	end
	
function h = imhist(im)
% compute image histogram. Display with bar(h)
  im = double(rgb2gray(im));
  im = uint8(round(im*256/max(im(:))));
  for l=0:255; h(l+1)=sum(im(:)==l); end
  
function im = imlogscale(im)
% enhance the intensity contrast using a log-scale
  im=double(im);
  im=log(1+im);
  im=uint8(im*255/max(im(:)));
  
function [x1,y1, m1, im1, sx, sy, iter] = max_and_zero(im1, dx, dy, iter)
  % max_and_zero: search for a maximum intensity point and zero image around it
  %
  % input:
  %   im1: image (gray m*n)
  %   dx,dy: area to zero after search
  % output:
  %   x1,y1: location in image
  %   m1:    intensity
  if nargin <2, dx=0.3; end
  if nargin <3, dy=0.3; end
  if nargin <4, iter=1; end
  if isscalar(dx) && dx<1
    dx = round(size(im1,1)*dx);
  end
  if isscalar(dy) && dy<1
    dy = round(size(im1,2)*dy);
  end
  
  % search for max intensity
  [m1,x1]=max(im1(:)); m1 = double(m1);
  [x1,y1]=ind2sub(size(im1),x1);
  [sx,sy, f1, f2] = peakwidth(im1, x1,y1);
  % blank image around max
  dx1 = round((x1-dx):(x1+dx));
  dy1 = round((y1-dy):(y1+dy));
  dx1=dx1(dx1>=1 & dx1 <=size(im1,1));
  dy1=dy1(dy1>=1 & dy1 <=size(im1,2));
  im1(dx1,dy1) = 0;
  
  if iter>20, return; end
  
  % recursive call if that guess is not acceptable
  % remove dead pixels (too sharp peaks) and image edges.
  if sx*sy < 4 || x1<5 || x1>size(im1,1)-4 || y1<5 || y1>size(im1,2)-4
    [x1,y1, m1, im1, sx, sy, iter] = max_and_zero(im1, dx, dy, iter+1);
    return
  end
  
  
function [s1,s2, f1, f2] = peakwidth(im, x,y, dx, dy)
% peakwidth: determine a peak width along X and Y
%
% input:
%   im:     image (rgb)
%   x,y:    peak location
%   dx,dy:  pixel size to extract around peak (optional, default=+/-5)

  if nargin < 4, dx=4; end
  if nargin < 5, dy=4; end

  % fist extract the portion of the image to analyze
  X=(x-dx):(x+dx); X=X(X>0 & X<size(im,1));
  Y=(y-dy):(y+dy); Y=Y(Y>0 & Y<size(im,2));
  
  % then project the image along dimension 1
  [s1, f1] = width1(Y, im(x,Y));
  [s2, f2] = width1(X, im(X,y));  

function [h,colors]=plot_im_points(im, t, x,y, m)
% plot_im_points: plot an image and the 'max' star locations
%
% input:
%   im:  image
%   t:   title (filename)
%   x,y: locations (vectors)
%   m:   intensity (vector)
  h = [];
  image(uint8(im)); title(t); hold on;
  colors =repmat('rgbmcy',1, numel(x));
  for p=1:numel(x)
    if ~x(p) || ~y(p) || m(p) <= 0 || ~isfinite(m(p)), continue; end
    h1=plot(y(p),x(p),colors(p));
    set(h1,'Marker','o', 'MarkerSize',10*m(p)/256);
    h = [ h h1 ];
  end
  hold off
  drawnow
  set(gcf, 'Name', t);

function plot_metrics(output)

  % plot some metrics from the transformation

  x=output.T(:,1); y=output.T(:,2);

  angle=output.angle;

  if numel(angle) == numel(output.time)
    t = cumsum(output.time);
    taxis = 'time [s]';
  else
    t = 1:numel(angle);
    taxis = 'image';
  end
  subplot(2,2,1);
  plot(x,y,'o:'); title('Translation. Start at (0,0)'); 
  xlabel('X [pix]'); ylabel('Y [pix]'); 

  subplot(2,2,2);
  plot(t,x,'ro', t,y,'bs'); title('Translation over time');
  xlabel(taxis); ylabel('X,Y/first image [pix]'); 

  try % this may fail if some images have different control point nb
    sx=cell2mat(cellfun(@(c)mean(c.sx), output.points,'UniformOutput',false)');
    sy=cell2mat(cellfun(@(c)mean(c.sy), output.points,'UniformOutput',false)');
    subplot(2,2,3);
    plot(t(2:end),sx.*sy); title('Control points area');
    xlabel(taxis); ylabel('Width^2 [pix^2]');
  end

  subplot(2,2,4);
  plot(t,angle); title('Rotation');
  xlabel(taxis); ylabel('Angle/first image [deg]'); 
  
function a = resolvefiles(a, t)
% resolvefiles: determine the fully qualified path of given names.
%   also resolves the wildcards expansion.
%
% output:
%   cellstr of filenames
  if nargin <2, t='Pick file(s)'; end
  if isempty(a)
    formats = imformats;
    
    % assemble the list of formats for the file selector
    f = {};
    f{end+1} = '*.*';
    f{end+1} = 'All Files (*.*)';
    f{end+1} = '*.jpg;*.JPG;*.jpeg;*.JPEG';
    f{end+1} = 'JPEG images (*.jpg)';
    f{end+1} = '*.png;*.PNG';
    f{end+1} = 'PNG images (*.jpg)';
    f{end+1} = '*.tif;*.TIF;*.tiff;*.TIFF';
    f{end+1} = 'TIFF images (*.tiff)';
    for index=1:numel(formats)
      this = formats(index);
      F=upper(this.ext);
      f{end+1} = [ sprintf('*.%s;', this.ext{:}) sprintf('*.%s;', F{:}) ];
      f{end+1} = this.description;
    end
    
    f = reshape(f, [ 2 index+4])';
   
    % get the files
    [a, pathname] = uigetfile( ...
          f, ...
          t, ...
          'MultiSelect', 'on');
    if isequal(a,0) || isempty(a), return; end
    a = strcat(pathname, a);
  else

    if ischar(a),    a = cellstr(a); end
    if iscell(a)
      new_a = {};
      for index=1:numel(a)
        this = a{index};
        if ischar(this)
          p    = fileparts(this);
          if isdir(this), p = this; end
          if isempty(p), p=pwd; end
          b    = dir(this); % in case we have a wildcard
          for n=1:numel(b)
            if ~b(n).isdir
              new_a{end+1} = fullfile(p, b(n).name);
            end
          end
        elseif isnumeric(this)
          new_a{end+1} = this;
        end
      end
      a = new_a;
    end
    
  end
  
function gray = rgb2gray (rgb)
    if ndims(rgb) == 2, gray = rgb; return; end
    gray = (0.2989 * rgb(:,:,1)) + (0.5870 * rgb(:,:,2)) + (0.1140 * rgb(:,:,3));
    gray = cast(gray, class(rgb));
    
function [R,t] = rigid_transform_3D(A, B)
  % http://nghiaho.com/?page_id=671
  % apply affine transform to coordinates with:
  %   B = (ret_R*A') + repmat(ret_t, 1, size(A,1))
    if nargin ~= 2
	    error('Missing parameters');
    end

    if ~(all(size(A) == size(B)))
      R=[]; t=[];
      return;
    end

    centroid_A = mean(A);
    centroid_B = mean(B);

    N = size(A,1);

    H = (A - repmat(centroid_A, N, 1))' * (B - repmat(centroid_B, N, 1));

    [U,S,V] = svd(H);

    R = V*U';

    if det(R) < 0 && size(V,2) > 3
        disp('Reflection detected');
        V(:,3) = -V(:,3); 
        R = V*U';
    end

    t = -R*centroid_A' + centroid_B';
    
function [s,f] = width1(x,s)
% width1: compute the width and central axis value of a 1D signal
%
% input:
%   x: axis (vector)
%   s: signal (vector)
%
% output:
%   s: second moment (width)
%   f: first moment (mean axis value)
  x=x(:); s=double(s(:));
  sum_s = sum(s);
  
  % first moment (mean)
  f = sum(s.*x)/sum_s; % mean value
  
  % second moment: sqrt(sum(x^2*s)/sum(s)-fmon_x*fmon_x);
  s = sqrt(sum(x.*x.*s)/sum_s - f*f);

function filename = write_stacked(im, filename, info)
  % we create a cell with PNG options for imwrite
    % clean the EXIF
    f = fieldnames(info);
    for index=1:numel(f)
      if isfield(info, f{index})
        this = info.(f{index});
        flag = ~isempty(this) && ~isstruct(this) && (ischar(this) || (isnumeric(this) && isscalar(this)));
        if ~flag
          info = rmfield(info, f{index});
        end
      end
    end
    % fix DateTime
    if isfield(info,'DateTime')
      try
        d = datenum(info.DateTime);
      catch
        index = find(info.DateTime == ':');
        if numel(index) == 4
          info.DateTime(index(1:2)) = '-';
        end
      end
    else info.DateTime = datestr(now);
    end
    desc = '';
    try
      desc = evalc('disp(info)');
    end
    if isempty(desc)
      try
        desc = disp(info);
      end
    end
    if ~isfield(info,'Software'),    info.Software    = mfilename; end
    if ~isfield(info,'XResolution'), info.XResolution = 350; end
    if ~isfield(info,'YResolution'), info.YResolution = 350; end
    if ~isfield(info,'Model'),       info.Model       = version; end
    if ~isfield(info,'Make'),        info.Make       = 'Matlab'; end
    args={ ...
      'XResolution', info.XResolution , ...
      'YResolution', info.YResolution , ...
      'Title', [ mfilename ' stacked images ']  , ...
      'Author', 'E. Farhi' , ...
      'Description', desc , ...
      'CreationTime', info.DateTime , ...
      'Software', info.Software , ...
      'Source', [ info.Model ' ' info.Make ] , ...
      'Comment', desc };
    % writing stacked image
    [p,f] = fileparts(filename);
    filename = fullfile(p, [ f '_stacked.png' ]);
    disp([ mfilename ': writing ' filename ])
    warning('off','MATLAB:writepng:changedTextChunk')
    try
      imwrite(im, filename, args{:});
    catch
      imwrite(im, filename);
    end
  


