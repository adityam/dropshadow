if not modules then modules = { } end modules ['t-dropshadow'] = {
    version   = 1.000,
    comment   = "Drop shadows",
    author    = "Peter Rolf",
    copyright = "Peter Rolf",
    email     = "peter.rolf@arcor.de",
    license   = "GNU General Public License"
}

local version_string = "version: 2016.01.18" -- for batch file

local lpeg = require("lpeg")
local lpegmatch = lpeg.match
local texdimen = tex.dimen
local tonumber = tonumber
local metafun = context.metafun
local showmessage = interfaces.showmessage
local unpack = table.unpack or unpack
local insert, fastcopy = table.insert, table.fastcopy
local format, todimen, lower, upper = string.format, string.todimen, string.lower, string.upper
local abs, floor, max, min, pow, sqrt = math.abs, math.floor, math.max, math.min, math.pow, math.sqrt
local rad, cos, sin = math.rad, math.cos, math.sin


local trace_drops = false  trackers.register("modules.dropshadow", function(v) trace_drops = v end)
-- tracker for path calculation (talkative!)
local trace_drops_path = false trackers.register("modules.dropshadow.path", function(v) trace_drops_path = v end)
local report_drops = logs.new("dropshadow")

thirddata       = thirddata       or { }
thirddata.drops = thirddata.drops or { }

local drops = thirddata.drops

drops.parameters = drops.parameters or { }
-- fallback parameter
local fallback = {
               background = "color",-- \framed parameter
          backgroundcolor = "white",-- \framed parameter
               colorspace = "",-- derived from shadowcolor, overriding is possible
                direction = "-45",
               fileformat = "png",
                   mppath = "",
                   offset = "1.69mm", -- 10px@150ppi
                pdistance = "0.51mm", -- 3px@150ppi
                 penumbra = "40",
                   psigma = "1.52mm", -- 9px@150ppi
                   radius = "0mm", -- boxshadow only
               resolution = "150",-- in ppi (pixel per inch)
                 rotation = "0",
                    setup = "",-- no setup in the defaults; never
    shadowbackgroundcolor = "white",
              shadowcolor = "black",
                udistance = "0mm",
                    umbra = "50",
                   usigma = "0.85mm", -- 5px@150ppi
}

-- also used as initial value for the defaults; these values can be changed by the user
drops.parameters.defaults = fastcopy(fallback)

-- mainly for the \dropsshowtable macro
-- must be initialised properly to handle too early calls in combination with \dropsshowtable...
-- note that these parameters are already interpreted (type is assigned)
-- the initialisation values don't count here, only 'type' must fit
drops.parameters.current  = {
          backgroundcolor = "white",
               colorspace = "",
                direction = -45,
               fileformat = "png",
                   mppath = "",
                   offset = "1.69mm",
                pdistance = "0.51mm",
                 penumbra = 40,
                   psigma = 9,
                   radius = 0,
               resolution = 150,
    shadowbackgroundcolor = "white",
              shadowcolor = "black",
                   shiftx = 0,
                   shifty = 0,
                    umbra = 50,
                udistance = "0mm",
                   usigma = 5,
                  xoffset = 7.2,
                  yoffset = -7.2,
}

-- http://lua-users.org/wiki/SimpleRound
local function round(num, idp)
    local mult = 10^(idp or 0)
    return floor(num * mult + 0.5) / mult
end

function drops.numberofpixels(n,idp)
    local pxdimen = tex.pdfpxdimen or tex.pxdimen
    return round(todimen(n)/pxdimen,idp)
end

local numberofpixels = drops.numberofpixels


function drops.resetdrops()
   context.setupdrops(fallback) -- a 'delayed' call even works when \setupdrops is still undefined
end

function drops.getdefault(name)
    local t = drops.parameters.defaults
    if t[name] then return t[name] end
end

function drops.getcurrent(name)
    local t = drops.parameters.current
    if t[name] then return t[name] end
end



-- OS specific
local im_batchfile, batchstart
local percent_sign, dir_separator, continue_line, exclamation_mark, left_parenthesis, right_parenthesis

-- Sets the maximal length of a '-draw' path (used two times per shadow) in the IM command line.
-- depends on OS; a command line on Windows can contain either 2047 or 8191 characters.
-- The readability is also a problem with very long paths.
local IM_MAX_PATH_LIMIT = 750 -- TODO: initial value, needs testing

local double_quote    = [["]]

if os.type == "windows" then
         im_batchfile = "drops-" .. tex.jobname .. ".bat"
    if trace_drops then
        im_batchstart = ":: 'drops' batch file for ImageMagick (www.imagemagick.org)\n" ..
                ":: " .. version_string .. "\n"
    else
        im_batchstart = "@ECHO OFF\n:: 'drops' batch file for ImageMagick (www.imagemagick.org)\n" ..
                ":: " .. version_string .. "\n"
    end
         comment_sign = [[::]]
         percent_sign = [[%%%%]]
        dir_separator = [[\]]
        continue_line = "^\n "
     exclamation_mark = [[!]]
     left_parenthesis = [[ ( ]]
    right_parenthesis = [[ ) ]]
else
         im_batchfile = "drops-" .. tex.jobname
        im_batchstart = "#!/bin/sh\n#'drops' batch file for ImageMagick (www.imagemagick.org)\n" ..
                "# " .. version_string .. "\n "
         comment_sign = [[#]]
         percent_sign = [[%%]]
        dir_separator = [[/]]
        continue_line = "\\\n "
     exclamation_mark = [[\!]]
     left_parenthesis = [[ \( ]]
    right_parenthesis = [[ \) ]]
end

-- store data of current shadow
local current_shadow_id
function drops.currentshadowid() context(current_shadow_id) end

function drops.spec()
    local c = drops.parameters.current
    return {
        id        = c.shadow_id, name     = c.shadow_name,
        umbra     = c.umbra    , usigma   = c.usigma,   udistance = c.udistance,
        penumbra  = c.penumbra , psigma   = c.psigma,   pdistance = c.pdistance,
        offset    = c.offset   , xoffset  = c.xoffset , yoffset   = c.yoffset,
        direction = c.direction, rotation = c.rotation,
        resolution = c.resolution,
        shadowcolor = c.shadowcolor,
        shadowbackgroundcolor = c.shadowbackgroundcolor,
        colorspace = c.colorspace,
        setup = c.setup
    }
end

function drops.showdropstable()
    local c = drops.spec()
    local px = tex.pdfpxdimen or tex.pxdimen
    local NC,VL,HL,AR = context.NC,context.VL,context.HL,context.AR
    context("\\vbox\\bgroup\\switchtobodyfont[6pt,ss]\\vbox\\bgroup")
    context.starttable({"|r|r|r|r|"})

    NC()
    context(format("\\hbox\\bgroup\\switchtobodyfont[4pt]graphic (%d)\\egroup",c.resolution))
    VL()
    context("opacity")NC()
    context("sigma")NC()
    context("distance (px)")
    AR()HL()

    NC()
    context("umbra")VL()
    context(c.umbra .. "\\letterpercent")NC()
    context(c.usigma .. "px")NC()
    context(format("%s",c.udistance))
    AR()

    NC()
    context("penumbra")VL()
    context(c.penumbra .. "\\letterpercent")NC()
    context(c.psigma .. "px")NC()
    context(format("%s",c.pdistance))
    AR()

    NC()VL()context.use(3)AR()

    NC()
    context("\\hbox\\bgroup\\switchtobodyfont[4pt]location rel.\\egroup")VL()
    context("direction")NC()
    context("offset")NC()
    context("xy-offset (px)")
    AR()HL()

    NC()
    context("value")
    VL()
    context(format("%d°",c.direction))NC()
    context(format("%.1fpx",numberofpixels(c.offset,1)))NC()
    context(format("%+.1f%+.1f",c.xoffset/px,c.yoffset/px))
    AR()

    NC()context.use(4)AR()

    NC()
    context("\\hbox\\bgroup\\switchtobodyfont[4pt]color\\egroup")NC()
    context.use(3)context(format("cs:%s, %s",c.colorspace,c.shadowcolor))
    AR()

    NC()
    context("\\hbox\\bgroup\\switchtobodyfont[4pt]file\\egroup")NC()
    context.use(3)
    context(drops.currentshadowname())
    AR()

    context.stoptable()
    context("\\egroup\\egroup")
end

local shadowdir = ".drops"
local shadowsubdir = tex.jobname or "nodir"
local shadowfileprefix = "sh"

-- initialisation on every run
-- update checksum of hash file
job.variables.tobesaved.DropsHashChecksum =
    file.checksum(format("%s/drops-%s.lua",shadowdir,tex.jobname)) or "no checksum"
-- update checksum of batch file
job.variables.tobesaved.DropsBatchChecksum =
    file.checksum(format("%s/%s",shadowdir,im_batchfile)) or "no checksum"



drops.filecounter = drops.filecounter or 0
drops.shadows = { }
drops.IMpaths = { }
local batch, run_batch_flag = { }, false
local handled = { }

-- the directory names are interpreted as macros, if "\" is used as directory separator;
-- "/" also works on MS-OS, so no need to mess around with catcodes
local function currentshadowpath() return (shadowdir .. "/" .. shadowsubdir) end
local current_shadow_name = "" -- initialisation is needed for \dropsshowtable
function drops.currentshadowname() context(currentshadowpath() .. "/" .. current_shadow_name) end


local IM_NO_VERSIONNUMBER     = -1
local IM_NO_VERSIONNUMBER_TXT = "not available"
-- get version number of ImageMagick
local function get_im_version()
    local identify_command = "identify -version" -- "identify -list configure" is not usable on Windows
    local C, Cs, P, R, S, V = lpeg.C, lpeg.Cs, lpeg.P, lpeg.R, lpeg.S, lpeg.V

    local cleanup = Cs( ( S("\n\r") / " " +1 )^0 )
    local versionstring = lpegmatch(cleanup,os.resultof(identify_command))

    local number, point, minus = R('09'), P('.'), P('-')
    local numbers = number^1
    local versionnr = C(numbers) * point * C(numbers) * point * C(numbers) * minus * C(numbers)
    local pattern   = P{ "ImageMagick " * versionnr + 1 * V(1) }

    local version, versiontext
    if versionstring ~= "" then
        local mainversion, subversion, subsubversion, buildversion = lpegmatch(pattern,versionstring)
        -- must be a number for comparision
        version = tonumber(format("%02d%02d%02d%02d", mainversion,subversion,subsubversion,buildversion))
        versiontext = format("%d.%d.%d-%d", mainversion,subversion,subsubversion,buildversion)
    else -- can't get version info
        showmessage("drops","noversioninfo")
        version, versiontext = IM_NO_VERSIONNUMBER, IM_NO_VERSIONNUMBER_TXT
    end
    report_drops("ImageMagick version = " .. versiontext)
    return version, versiontext
end




local function ensure_dir(path,dir)
    if not io.exists(path) then -- io.exists() works only with an absolute path
        local state,errortxt = lfs.mkdir(path)
        if state ~= nil then
            if trace_drops then showmessage("drops","createdir",dir) end
        else if lower(errortxt) ~= "file exists" then
            showmessage("drops","direrror",{dir,errortxt}) end
        end
    end
end


local function load_hash()
    local absolute_shadowpath = lfs.currentdir() .. dir_separator .. shadowdir
    local hashfile = absolute_shadowpath .. dir_separator .. "drops-" .. tex.jobname .. ".lua"

    if io.exists(hashfile) then lua.registercode(hashfile)
    --    inspect(thirddata.drops.shadows)
    end
end

load_hash() -- at every run
-- get version number of IM; if this fails, everything else is in vain
if not drops.im_version then
    drops.im_version, drops.im_versiontext = get_im_version() end
if drops.im_version == IM_NO_VERSIONNUMBER then context("\\dropsoff") end -- no working IM, no drops


local function save_hash() -- at every \bye (see batch_control() )
    -- check directory structure
    local absolute_shadowpath = lfs.currentdir() .. dir_separator .. shadowdir
    ensure_dir(absolute_shadowpath,shadowdir)
    absolute_subpath = absolute_shadowpath .. dir_separator .. shadowsubdir
    ensure_dir(absolute_subpath,shadowsubdir)

    local hashfile = absolute_shadowpath .. dir_separator .. "drops-" .. tex.jobname .. ".lua"

    local imhash, err = io.open(hashfile,"w")

    if imhash then
        imhash:write(format("thirddata.drops.im_version=%d\n",drops.im_version or IM_NO_VERSIONNUMBER))
        imhash:write(format("thirddata.drops.im_versiontext=\"%s\"\n",drops.im_versiontext or "not available"))
        imhash:write(format("thirddata.drops.filecounter=%d\n",drops.filecounter or 0))
        imhash:write(table.serialize(drops.shadows,"thirddata.drops.shadows"))
        imhash:write("\n")
        imhash:write(table.serialize(drops.IMpaths,"thirddata.drops.IMpaths"))
        imhash:close()
        job.variables.tobesaved.DropsHashChecksum =
            file.checksum(format("%s/drops-%s.lua",shadowdir,tex.jobname)) or "no checksum"
    else
        showmessage("drops","ioopenerror",{hashfile,err})
    end
end


local function batch_start()
    -- check directory structure
    local absolute_shadowpath = lfs.currentdir() .. dir_separator .. shadowdir
    ensure_dir(absolute_shadowpath,shadowdir)
    absolute_subpath = absolute_shadowpath .. dir_separator .. shadowsubdir
    ensure_dir(absolute_subpath,shadowsubdir)

    local batchfile = absolute_shadowpath .. dir_separator .. im_batchfile

    local imbatch, err = io.open(batchfile,"w")

    if imbatch then
        if im_batchstart ~= "" then
            imbatch:write(im_batchstart,"\n")
        end
        for i,cmd in pairs(batch) do
            imbatch:write(cmd,"\n")
        end
        imbatch:close()
        showmessage("drops","startbatch",batchfile)
        os.execute("sh " .. batchfile)
    else
        showmessage("drops","ioopenerror",{batchfile,err})
    end
    job.variables.tobesaved.DropsBatchChecksum =
        file.checksum(format("%s/%s",shadowdir,im_batchfile)) or "no checksum"
end

function drops.batch_control()
    save_hash() -- save new hash; a different checksum forces a new run, so that the created graphics can be included
    if run_batch_flag then batch_start() end
end



local function eight_bit(n)
    return(round(n*255))
end

function drops.color2text(name)
    local ctxt, model
    if type(attributes.colors.spec) == "function" then -- introduced in beta 04.09.2012
        local t  = attributes.colors.spec(name)
        model = t.model

        local eb = eight_bit
        if model == "gray" then
            ctxt = format("%s(%s)",model,eb(t.s))
        elseif model == "rgb" then
            ctxt = format("%s(%s,%s,%s)",model,eb(t.r),eb(t.g),eb(t.b))
        elseif model == "cmyk" then
            ctxt = format("%s(%s,%s,%s,%s)",model,eb(t.c),eb(t.m),eb(t.y),eb(t.k))
        elseif model == "spot" then -- ???? untested
            ctxt = format("%s(%s,%s,%s,%s)",model,eb(t.c),eb(t.m),eb(t.y),eb(t.k))
            model = "cmyk"
        else ctxt = name ; model = "rgb" -- just in case
        end
    else -- use X11 colors otherwise (ConTeXt colors are not supported in older versions)
        ctxt = name ; model = "rgb"
    end
    return ctxt,model
end



-- no need to force the usage of full pixel for the offset
function drops.locateshadow(angle,offset,rotation,shiftx,shifty)
    local a,o,sx,sy,x,y
    a = rad(tonumber(angle) - tonumber(rotation))
    o = todimen(offset)
    sx = tonumber(shiftx)*65536
    sy = tonumber(shifty)*65536
    x,y = cos(a)*o+sx, sin(a)*o-sy
    texdimen.dropsXPos, texdimen.dropsYPos = x,-y -- write result into register
end


local TIME_INTERVAL = 1/25

local PATH_MIN_OFFSET = 1 -- minimal offset (in px)
local PATH_MAX_ERROR = 1.5 -- maximal allowed difference before complaining (sqrt(2) for a rectangular corner; spiky ones have much higher values)

local FACTOR_POST = 1/3
local FACTOR_PRE  = 2/3

local ACCURACY_POINT = 5 -- accuracy for points
local ACCURACY_SHIFT = 3 -- correction shift
local ACCURACY_FX = ACCURACY_POINT -- accuracy for the function parameters m and n

-- TODO: these values base on very few examples; needs a lot more testing for robust values

-- threshold to detect collinearity
local AREA_THRESHOLD = 0.00001
-- threshold to detect parallelism
local PARALLELISM_THRESHOLD = 0.00050

-- some helper macros for point locations in the path
local function point(p,i)  return { p[i][1], p[i][2] } end
local function cp1(p,i)    return { p[i][3], p[i][4] } end -- pre control point
local function cp2(p,i)    return { p[i][5], p[i][6] } end -- post


local function calculate_point(A,B,C,D,t)
    local ax,ay,bx,by,cx,cy,dx,dy = A[1],A[2],B[1],B[2],C[1],C[2],D[1],D[2]

    if     t == 0 then return ax,ay
    elseif t == 1 then return dx,dy
    else
        local t2,t3,tx,ty
        t2 = t*t
        t3 = t*t2
        -- http://en.wikipedia.org/wiki/B%C3%A9zier_curve
        -- $ C(t) = (-P_0+3P_1-3P_2+P_3)t^3 + (3P_0-6P_1+3P_2)t^2 + (-3P_0+3P_1)t + P_0, t\elem[0,1] $
        -- this formula is taken from the German Wikipedia (not part of the English version)
        -- http://de.wikipedia.org/wiki/B%C3%A9zierkurve#Kubische_B.C3.A9zierkurven_.28n.3D3.29
        -- tx = (-ax+3*bx-3*cx+dx)*t3 + (3*ax-6*bx+3*cx)*t2 + (-3*ax+3*bx)*t + ax
        -- ty = (-ay+3*by-3*cy+dy)*t3 + (3*ay-6*by+3*cy)*t2 + (-3*ay+3*by)*t + ay
        tx = (-ax+3*(bx-cx)+dx)*t3 + (3*(ax-2*bx+cx))*t2 + (3*(-ax+bx))*t + ax -- saves 3 multiplications (*2)
        ty = (-ay+3*(by-cy)+dy)*t3 + (3*(ay-2*by+cy))*t2 + (3*(-ay+by))*t + ay

        return tx,ty
    end
end

local function get_local_error(A,B,C,D,E,F,G,H,d,i)
    local interval = i or TIME_INTERVAL
    local maxerror,locerror,time = 0
    local px,py,qx,qy

    for t=0,1,interval do -- TODO: what values are needed for acceptable results? equidistant or a list of fixed values?
        px,py = calculate_point(A,B,C,D,t)
        qx,qy = calculate_point(E,F,G,H,t)
        locerror = sqrt(pow(qy-py,2) + pow(qx-px,2)) -abs(d)

        if abs(locerror) > maxerror then maxerror = locerror; time = t end
    end

    return maxerror,time
end


-- support for passvariable() and passarrayvariable()
local function is_path_array(pd)
    return (type(pd[1][1]) == "table")
end

-- and multi paths
local function is_multi_path(pd)
    return (is_path_array(pd) and (#pd>1))
end

-- relies on multi path data structure
local function is_closed_path(pd)
    for i = 1,#pd do
        local n = #pd[i]
        if not (pd[i][1][1] == pd[i][n][1] and pd[i][1][2] == pd[i][n][2]) then
            return false
        end
    end
    return true
end

-- based on the ** CODE 1 ** example at
-- http://stackoverflow.com/questions/2587751/an-algorithm-to-find-bounding-box-of-closed-bezier-curves
local function get_single_path_extrema(path)
    local p = fastcopy(path); local n = #p
    local P1,P2,P3,P4 = point(p,1)
    local r = { {P1[1]} , {P1[2]} } -- store possible extrema (roots) for x and y; start point is first candidate
    local rx,ry = r[1], r[2] -- separate arrays for the x|y coordinates

    local a,b,c,d,sqrtd,t
    local function f(t,i) -- same formula as in calculate_point(), but only for one coordinate axis
        return (-P1[i]+3*(P2[i]-P3[i])+P4[i])*pow(t,3) + (3*(P1[i]-2*P2[i]+P3[i]))*pow(t,2) + (3*(-P1[i]+P2[i]))*t + P1[i]
    end

    local ZERO_THRESHOLD = pow(10,-12) -- a more 'stable' test for zero; coordinates aren't that accurate anyway...
    local i,k

    for i = 1,n-1 do
        P2,P3,P4 = cp2(p,i), cp1(p,i+1), point(p,i+1)

        insert(rx,P4[1]); insert(ry,P4[2]) -- start|end points are always candidates

        for k=1,2 do
            a = -3*P1[k]  +9*P2[k] -9*P3[k] +3*P4[k]
            b =  6*P1[k] -12*P2[k] +6*P3[k]
            c = -3*P1[k]  +3*P2[k]

            if abs(a) < ZERO_THRESHOLD then -- a == 0
                if not (abs(b) < ZERO_THRESHOLD) then -- not (b == 0)
                    t = -c / b
                    if ( 0 < t and t < 1) then insert(r[k],f(t,k)) end
                end
            else
                d = pow(b,2) - 4*c*a ; sqrtd = sqrt(d)
                if not (d < 0) then
                    t = (-b + sqrtd)/(2*a)
                    if ( 0 < t and t < 1 ) then insert(r[k],f(t,k)) end
                    t = (-b - sqrtd)/(2*a)
                    if ( 0 < t and t < 1) then insert(r[k],f(t,k)) end
                end
            end
        end -- x,y

        P1 = P4
     end -- path

     local x_min,x_max,y_min,y_max
     x_min = min(unpack(rx)) ; x_max = max(unpack(rx))
     y_min = min(unpack(ry)) ; y_max = max(unpack(ry))

     if trace_drops_path then
         report_drops("get_single_path_extrema:  x=(%.5f,%.5f), y=(%.5f,%.5f)",x_min,x_max,y_min,y_max)
     end
     return x_min,x_max,y_min,y_max
end


function drops.get_path_extrema(p)
    local xmin,ymin = p[1][1][1], p[1][1][2] -- point 1 of path 1
    local xmax,ymax = xmin,ymin
    local txmin,txmax,tymin,tymax
    for i=1,#p do -- get extrema for every sub path
        txmin,txmax,tymin,tymax = get_single_path_extrema(p[i])
        if txmin < xmin then xmin = txmin end
        if txmax > xmax then xmax = txmax end
        if tymin < ymin then ymin = tymin end
        if tymax > ymax then ymax = tymax end
    end

    return xmin,xmax,ymin,ymax
end


function drops.get_path_minima(path)
    local xmin,xmax,ymin,ymax = drops.get_path_extrema(path)
    return xmin,ymin
end

function drops.get_path_maxima(path)
    local xmin,xmax,ymin,ymax = drops.get_path_extrema(path)
    return xmax,ymax
end




-- http://mathworld.wolfram.com/Collinear.html
local function points_are_collinear(A,B,C,D)
    local ax,ay,bx,by,cx,cy,dx,dy = A[1],A[2],B[1],B[2],C[1],C[2],D[1],D[2]

    local t1,t2
    -- the area of a triangle is sero, if all points are collinear
    t1 = ax*(by-cy)+bx*(cy-ay)+cx*(ay-by)
    t2 = bx*(cy-dy)+cx*(dy-by)+dx*(by-cy)

    if abs(t1) < AREA_THRESHOLD and abs(t2) < AREA_THRESHOLD then
        return true
    else
        return false
    end
end



-- http://pomax.github.io/bezierinfo/#intersections
local function get_isp(A,B,C,D)
    local ax,ay,bx,by,cx,cy,dx,dy = A[1],A[2],B[1],B[2],C[1],C[2],D[1],D[2]
    local isp = nil

    local d = (ax-bx)*(cy-dy)-(ay-by)*(cx-dx)
--    if d == 0 then
    -- parallelism test (not that robust; very small values can result in very high values for the isp coordinates)
    if abs(d) > PARALLELISM_THRESHOLD then
        local t1,t2,nx,ny = ax*by-ay*bx, cx*dy-cy*dx
        nx = t1*(cx-dx)-(ax-bx)*t2
        ny = t1*(cy-dy)-(ay-by)*t2
        isp = { nx/d, ny/d }
    end

    if trace_drops_path then
        if isp then
            report_drops("get_isp: returning (%.1f,%.1f)",isp[1],isp[2])
        else
            report_drops("get_isp: returning nil")
        end
    end

    return isp
end



-- some helper functions
local function line_length(A,B)
    return sqrt(pow(B[2]-A[2],2)+pow(B[1]-A[1],2)) -- use the Pythagorean theorem (gradient is a right-angled triangle)
end

local function unit_vector(A,B,length)
    local l = length or line_length(A,B)
    if round(l,ACCURACY_POINT) > 0 then
        return { -(B[2]-A[2])/l, (B[1]-A[1])/l }
    else
        return { 0,0 }
    end
end

-- algorithm by Tiller-Hanson;  based on the pseudocode written in
-- https://groups.google.com/forum/?_escaped_fragment_=searchin/comp.graphics.algorithms/#!search/outsetting$20insetting$20a$20bezier$20path/comp.graphics.algorithms/Ec6lEX_ppHI/do6RoQsEbfYJ
-- calculates both offsets (umbra,penumbra) to save some runtime
local function get_offsets(A,B,C,D,d1,d2)
    local d1_distance, d2_distance = abs(d1) >= PATH_MIN_OFFSET, abs(d2) >= PATH_MIN_OFFSET

    if d1_distance or d2_distance then
        local x,y = 1,2 -- for readability
        local ax,ay,bx,by,cx,cy,dx,dy = A[x],A[y],B[x],B[y],C[x],C[y],D[x],D[y]

        local length_AB,length_BC,length_CD -- line length (needed for translation into unit vectors)
        length_AB = line_length(A,B)
        length_BC = line_length(B,C)
        length_CD = line_length(C,D)

        if trace_drops_path then
            report_drops("get_offsets: (A,B,C,D),d1,d2 = ( (%d,%d), (%d,%d), (%d,%d), (%d,%d) ),%d,%d",
                     ax,ay,bx,by,cx,cy,dx,dy,d1,d2)
            report_drops("get_offsets: length_AB=%.1f, length_BC=%.1f, length_CD=%.1f",
                     length_AB,length_BC,length_CD)
        end

        local AB,BC,CD = length_AB>0, length_BC>0, length_CD>0

        -- get normals for the three sections and divide them by their length; we get unit vectors
        local unit_AB,unit_BC,unit_CD
        -- x=-y, y=x
        if     AB then unit_AB = unit_vector(A,B,length_AB) -- A~=B
        elseif BC then unit_BC = unit_vector(B,C,length_BC) -- A=B, B~=C
            unit_AB = unit_BC
        elseif CD then unit_CD = unit_vector(C,D,length_CD) -- A=B=C, C~=D
            unit_AB = unit_CD
        else -- all four points are equal (case is finally handled here, no further check)
            if trace_drops_path then
                report_drops("get_offsets: *warning: all four points are equal; returning the original path")
            end
            return A,A,A,A,A,A,A,A -- no calculation possible
        end

        if not unit_BC then
            if     BC then unit_BC = unit_vector(B,C,length_BC) -- B~=C
            elseif AB then unit_BC = unit_AB -- B=C, A~=B
            else --if CD then
                unit_CD = unit_vector(C,D,length_CD) -- A=B=C, C~=D
                unit_BC = unit_CD
            end
        end

        if not unit_CD then
            if     CD then unit_CD = unit_vector(C,D,length_CD) -- C~=D
            elseif BC then unit_CD = unit_BC -- C=D, B~=C
            else --if AB then
                unit_CD = unit_AB -- B=C=D, A~=B
            end
        end
        -- all the calculations with the original path are done now...

        -- temporal points for the shifted lines (A,B)->(T1,T2),  (B,C)->(T3,T4), (C,D)->(T5,T6)
        local T1,T2,T3,T4,T5,T6
        local tx,ty

        local P1,P2,P3,P4 -- points for the umbra offset
        if d1_distance then
            -- now shift the original lines into the direction of the unit vectors
            tx,ty = unit_AB[x]*d1, unit_AB[y]*d1; T1 = { ax+tx, ay+ty }; T2 = { bx+tx, by+ty }
            tx,ty = unit_BC[x]*d1, unit_BC[y]*d1; T3 = { bx+tx, by+ty }; T4 = { cx+tx, cy+ty }
            tx,ty = unit_CD[x]*d1, unit_CD[y]*d1; T5 = { cx+tx, cy+ty }; T6 = { dx+tx, dy+ty }

            P1,P2,P3,P4 = T1,{},{},T6 -- first and last points are ok (for now);

            if points_are_collinear(T1,T2,T3,T4) then

                tx = T6[x]-T1[x]; ty = T6[y]-T1[y]
                P2 = { T1[x]+ FACTOR_POST*tx, T1[y]+ FACTOR_POST*ty }
                P3 = { T1[x]+ FACTOR_PRE *tx, T1[y]+ FACTOR_PRE *ty }

                if trace_drops_path then
                    report_drops("get_offsets: 1-4 are collinear; setting P2(post) to (%.1f,%.1f), P3(pre) to (%.1f,%.1f)",P2[x],P2[y],P3[x],P3[y])
                end
            else
                P2 = get_isp(T1,T2,T3,T4)
                if not(P2) then
                    P2 = { T1[x] + FACTOR_POST*(T4[x]-T1[x]), T1[y] + FACTOR_POST*(T4[y]-T1[y]) }
                    if trace_drops_path then
                        report_drops("get_offsets: 1-3 are collinear; setting P2(post) to (%.1f,%.1f)",P2[x],P2[y])
                    end
                end
                P3 = get_isp(T3,T4,T5,T6)
                if not(P3) then
                    P3 = { T3[x] + FACTOR_PRE*(T6[x]-T3[x]), T3[y] + FACTOR_PRE*(T6[y]-T3[y]) }
                    if trace_drops_path then
                        report_drops("get_offsets: 2-4 are collinear; setting P3(pre) to (%.1f,%.1f)",P3[x],P3[y])
                    end
                end

            end

        else P1,P2,P3,P4 = A,B,C,D
        end



        local Q1,Q2,Q3,Q4 -- points for the penumbra offset
        if d1==d2 then
            Q1,Q2,Q3,Q4 = P1,P2,P3,P4
        elseif d2_distance then
            -- now shift the original lines into the direction of the unit vectors
            tx,ty = unit_AB[x]*d2, unit_AB[y]*d2; T1 = { ax+tx, ay+ty }; T2 = { bx+tx, by+ty }
            tx,ty = unit_BC[x]*d2, unit_BC[y]*d2; T3 = { bx+tx, by+ty }; T4 = { cx+tx, cy+ty }
            tx,ty = unit_CD[x]*d2, unit_CD[y]*d2; T5 = { cx+tx, cy+ty }; T6 = { dx+tx, dy+ty }

            Q1,Q2,Q3,Q4 = T1,{},{},T6 -- first and last points are ok (for now);

            if points_are_collinear(T1,T2,T3,T4) then

                tx = T6[x]-T1[x]; ty = T6[y]-T1[y]
                Q2 = { T1[x]+ FACTOR_POST*tx, T1[y]+ FACTOR_POST*ty }
                Q3 = { T1[x]+ FACTOR_PRE *tx, T1[y]+ FACTOR_PRE *ty }

                if trace_drops_path then
                    report_drops("get_offsets: 1-4 are collinear; setting Q2(post) to (%.1f,%.1f), Q3(pre) to (%.1f,%.1f)",Q2[x],Q2[y],Q3[x],Q3[y])
                end
            else
                Q2 = get_isp(T1,T2,T3,T4)
                if not(Q2) then
                    Q2 = { T1[x] + FACTOR_POST*(T4[x]-T1[x]), T1[y] + FACTOR_POST*(T4[y]-T1[y]) }
                    if trace_drops_path then
                        report_drops("get_offsets: 1-3 are collinear; setting Q2(post) to (%.1f,%.1f)",Q2[x],Q2[y])
                    end
                end
                Q3 = get_isp(T3,T4,T5,T6)
                if not(Q3) then
                    Q3 = { T3[x] + FACTOR_PRE*(T6[x]-T3[x]), T3[y] + FACTOR_PRE*(T6[y]-T3[y]) }
                    if trace_drops_path then
                        report_drops("get_offsets: 2-4 are collinear; setting Q3(pre) to (%.1f,%.1f)",Q3[x],Q3[y])
                    end
                end

            end

        else Q1,Q2,Q3,Q4 = A,B,C,D
        end

        if trace_drops_path then
            report_drops("get_offsets: returning ( (%.1f,%.1f), (%.1f,%.1f), (%.1f,%.1f), (%.1f,%.1f) )",
                         P1[x],P1[y],P2[x],P2[y],P3[x],P3[y],P4[x],P4[y])
            report_drops("get_offsets:           ( (%.1f,%.1f), (%.1f,%.1f), (%.1f,%.1f), (%.1f,%.1f) )",
                         Q1[x],Q1[y],Q2[x],Q2[y],Q3[x],Q3[y],Q4[x],Q4[y])
        end


        return P1,P2,P3,P4,Q1,Q2,Q3,Q4


    else
        if trace_drops_path then
            report_drops("get_offsets: d1=d2=0; returning original points")
        end
        return A,B,C,D,A,B,C,D
    end

end


-- start mask graphic at (x,y) and not necessarily at (0,0); one free pixel on each side is preferable;
-- this avoids ugly looking border cases (alpha mask is visible) when rendering
local default_x_offset,default_y_offset = 1,1

-- 'normalize' path by shifting it;
-- all path points with negative coordinates are not part of the IM canvas (in other words: they are not visible)!
-- control points can be negative though (as long as they are not lying on the path)
local function path_shift_to(p,xoffset,yoffset,xmin,ymin)
    local x_offset,y_offset = xoffset or default_x_offset, yoffset or default_y_offset
    local x_min,y_min,sx,sy
    if xmin and ymin then
        x_min,y_min = xmin,ymin
    else
        x_min,y_min = drops.get_path_minima(p)
    end

    sx = -x_min +x_offset; sy = -y_min +y_offset -- the final shift values for the path
    if trace_drops_path then
        report_drops("path_sto: minima: (%.5f, %.5f), default offset: (%d,%d)",x_min,y_min,x_offset,y_offset)
        report_drops("               shifting path by (%.5f, %.5f)",sx,sy)
    end

    local i,j,k
    for i=1,#p do -- for every sub path
        q = p[i]
        -- substract the minima from every coordinate (control points may get negative coordinates)
        for j=1,#q do
            for k=0,2 do
                q[j][(2*k)+1] = round(q[j][(2*k)+1]+sx,ACCURACY_POINT)
                q[j][(2*k)+2] = round(q[j][(2*k)+2]+sy,ACCURACY_POINT)
            end
        end
    end

    return sx,sy -- return the shift values
end



-- for the default boxshadow
function drops.getboxshadowid(width,height,radius)
    local w,h,r = numberofpixels(width), numberofpixels(height), numberofpixels(radius)
    local name = format("drops:boxshadow-w%dh%dr%d",w,h,r)
    return name
end


local function generate_default_path(id,width,height,radius)
    local w,h,r
    -- I directly use pixel numbers here to save the path scaling (bp to px) for IM.
    -- Normally you cannot assume, that a given MP path is already scaled correctly for IM usage.
    -- This is the only exception, as extra handling is needed in the drops.shadow() routine!
    -- If you use a path that is given in 'bp' units (the default case), you have to scale it for IM ('px' based, depends on resolution value).
    -- Use something like this in your path generation function:
    --     w = number.tobasepoints(width);  h = number.tobasepoints(height); ...

    w = numberofpixels(width); h = numberofpixels(height); r = numberofpixels(radius)

    metafun.start()
    metafun("save p,pid; path p; string pid;")
    metafun("pid := \"%s\";",id)
    --metafun("show pid;")
    if r == 0 then
        metafun("p:= unitsquare xyscaled(%s,%s);",w,h) -- five path segments
    elseif r >= 15 then -- TODO: add as user parameter and make some tests for a good default value
        metafun("primarydef p wellrounded d = (hide(path __helpcircle__; __helpcircle__ := unitcircle scaled d;) subpath(4,6) of __helpcircle__ -- subpath(6,8) of __helpcircle__ xshifted(xpart(lrcorner p)-d) -- subpath(0,2) of __helpcircle__ shifted(xpart(lrcorner p)-d,ypart(urcorner p)-d) -- subpath(2,4) of __helpcircle__ yshifted(ypart(ulcorner p)-d) -- cycle) enddef ; p:= unitsquare xyscaled(%s,%s) wellrounded %s;",w,h,r) -- thirdteen path segments
    else
        metafun("p:= unitsquare xyscaled(%s,%s) smoothed %s;",w,h,r) -- nine path segments (bad results with bigger radii)
    end
    metafun("passvariable(pid,p);")
    metafun("setbounds currentpicture to boundingbox(p);") -- deal with the size (box it, if unwanted)
    metafun.stop()
end

function drops.createboxshadowpath(width,height,radius)
    local pathid = drops.getboxshadowid(width,height,radius)
    if not metapost.variables[pathid] then
        generate_default_path(pathid,width,height,radius)
--        context(function() inspect(metapost.variables[pathid]) end)
    end
end

-- scale, swap y-axis and shift to get positive only path coordinates
local function adapt_path_to_IM(name,pathid,resolution) -- name of original MP path, name of scaled IM path
    if metapost.variables[name] then
        local p = fastcopy(metapost.variables[name])

        -- scale from 'bp' to current drops resolution for already predefined MP paths only
        -- swap orientation of the y-axis (the resulting shift can be ignored; we clean up later)
        local s = resolution and resolution/72 or 1

        if not is_path_array(p) then p = { p } end -- use a common data structure

        if not is_closed_path(p) then
            showmessage("drops","noclosedpath",name)
            -- return nil
        end

        local i,j,q
        if s ~= 1 then
            for i=1,#p do -- for all sub paths
                q = p[i]
                for j=1,#q do
                    q[j][1] =  q[j][1] * s; q[j][3] =  q[j][3] * s; q[j][5] =  q[j][5] * s -- x-values
                    q[j][2] = -q[j][2] * s; q[j][4] = -q[j][4] * s; q[j][6] = -q[j][6] * s -- swap orientation of y-values
                end
            end
        else -- special handling of boxshadow paths (already in px)
            for i=1,#p do -- for all sub paths
                q = p[i]
                for j=1,#q do
                    q = p[i]
                    q[j][2] = -q[j][2]; q[j][4] = -q[j][4]; q[j][6] = -q[j][6] -- swap orientation of y-values
                end
            end
        end

        path_shift_to(p)

        -- path should be fine now (right scaled and only positive coordinates)
        -- store the 'IM ready' path for later usage
        drops.IMpaths[pathid] = {
            path = fastcopy(p)
        }
        return true
    else
        showmessage("drops","nomppath",name)
        return nil
    end
end

-- creates '-draw' paths for IM
local function get_IM_paths(name,udistance,pdistance)

    local ud,pd = udistance, pdistance

    if trace_drops_path then
        report_drops("get_IM_paths: ******************************************")
        report_drops("get_IM_paths: name=%s, udistance=%s, pdistance=%s",name,ud,pd)
    end

    local IMpaths = drops.IMpaths
    local resolution = drops.parameters.current.resolution
    local pathid = name .. "_" .. resolution .. "ppi"

    local impath = IMpaths[pathid] and true or adapt_path_to_IM(name,pathid,resolution)

    if impath then
        local x,y = 1,2
        local u,p = {},{} -- umbra/penumbra path array

        local A,B,C,D
        local a,apre,apost,b,bpre,bpost -- variables for u
        local c,cpre,cpost,d,dpre,dpost -- variables for p
        local i,j -- indices for start|end point
        local isp -- intersection point
        local error,errortime

        local sp -- counter for sub paths
        local nofsp = #IMpaths[pathid].path

        for sp=1,nofsp do
            local o = IMpaths[pathid].path[sp] ; local nofpoints = #o

            u[sp],p[sp] = {},{}
            local oldb,oldbpre, oldd,olddpre
            -- get 'pre' control value of the first point, so that we can store it in the loop with both control points
            a,apost,oldbpre,oldb, c,cpost,olddpre,oldd = get_offsets(point(o,nofpoints),cp2(o,nofpoints),cp1(o,1),point(o,1),ud,pd)

            for i = 1,nofpoints do
                if trace_drops_path then
                    report_drops("get_IM_paths:  ------------------------ %d [%d] ------------------------",i,sp)
                end

                j = i%nofpoints +1 -- needed for closed path
                A,B,C,D = point(o,i),cp2(o,i),cp1(o,j),point(o,j)
                a,apost,bpre,b, c,cpost,dpre,d = get_offsets(A,B,C,D,ud,pd)

                -- 'a' and 'oldb' have the same origin, but are used in a different path segment.
                -- 'oldb' is used as end point (+pre) and 'a' as starting point (+post)
                -- Currently the isp of oldbpre-oldb and a-apost is used as new point; the number of points is unchanged,
                -- no path intersection with negative offsets.
                -- Better, but much harder to implement: create an arc between the points; number of points is doubled
                -- and intersection problems.
                isp = get_isp(a,apost,oldbpre,oldb)

                if isp then
                    if trace_drops_path then
                        error,errortime = get_local_error(A,B,C,D,a,apost,bpre,b,ud)
                        if error > PATH_MAX_ERROR and i ~= nofpoints then -- fishy error values with last point; why?
                            report_drops("get_IM_paths: *warning: umbra max_error=%.5f (t=%.2f)",error,errortime)
                        end
                    end
                    insert(u[sp], { isp[x],isp[y], oldbpre[x],oldbpre[y], apost[x],apost[y] })
                else
                    insert(u[sp], { oldb[x],oldb[y], oldbpre[x],oldbpre[y], apost[x],apost[y] })
                    if trace_drops_path then
                        report_drops("get_IM_paths: failed to calculate ISP (umbra)")
                    end
                end
                isp =  get_isp(c,cpost,olddpre,oldd)
                if isp then
                    if trace_drops_path then
                        error,errortime = get_local_error(A,B,C,D,c,cpost,dpre,d,pd)
                        if error > PATH_MAX_ERROR and i ~= nofpoints then -- fishy error values with last point; why?
                            report_drops("get_IM_paths: *warning: penumbra max_error=%.5f (t=%.2f)",error,errortime)
                        end
                    end
                    insert(p[sp],{ isp[x],isp[y], olddpre[x],olddpre[y], cpost[x],cpost[y] } )
                else
                    insert(p[sp],{ oldd[x],oldd[y], olddpre[x],olddpre[y], cpost[x],cpost[y] } )
                    if trace_drops_path then
                        report_drops("get_IM_paths: failed to calculate ISP (penumbra)")
                    end
                end

                oldb,oldbpre,oldd,olddpre = b,bpre,d,dpre
            end

        end


        -- The center point of an offset path is not necessarily the center of the origin path.
        -- We must compensate that with a correction shift of the shadow bitmap.
        -- This is another drawback of the used method.
        local uxmin,uxmax,uymin,uymax = drops.get_path_extrema(u)
        local pxmin,pxmax,pymin,pymax
        if ud==pd then
            pxmin,pxmax,pymin,pymax = uxmin,uxmax,uymin,uymax
        else
            pxmin,pxmax,pymin,pymax = drops.get_path_extrema(p)
        end

        -- calculate the needed shift to "recenter" the graphic
        -- penumbra path should always be bigger than umbra path (fails otherwise)
        local divisor = 2*(resolution/72.72) -- shift only half of the difference and also consider the current resolution
        local shiftx = round((pxmin-uxmin + pxmax-uxmax)/divisor,ACCURACY_SHIFT)
        local shifty = round((pymin-uymin + pymax-uymax)/divisor,ACCURACY_SHIFT)
        if trace_drops_path then
            report_drops("get_IM_paths: shift = (%.5f,%.5f)",shiftx,shifty)
        end

        -- shift the bigger penumbra area to get positive only coordinates and capture the shift values
        local px,py = path_shift_to(p,default_x_offset,default_y_offset,pxmin,pymin)
        -- now we shift the smaller umbra shadow on the bigger penumbra canvas, so that they "match" again
        -- shift the umbra path the same way as the penumbra path and add its own correction terms
        path_shift_to(u,px+uxmin,py+uymin,uxmin,uymin)

        -- calculate the correct bounding box
        -- +2 (one free pixel on each side) +1 (canvas coordinates start at (1,1) )
        local bbx,bby = pxmax-pxmin+3, pymax-pymin+3

        -- now create the IM paths...
        local  upath, ppath, close_path = "", "", ""
        local uprev, pprev, ucurrent, pcurrent

        for sp=1,nofsp do
            uprev, pprev = nil, nil

            for i=1,#u[sp] do
                ucurrent = u[sp][i]
                if not uprev then
                    upath = upath .. string.format("M %d,%d",ucurrent[1],ucurrent[2])
                else
                    upath = upath .. string.format(" C %d,%d %d,%d %d,%d",
                        uprev[5],uprev[6],ucurrent[3],ucurrent[4],ucurrent[1],ucurrent[2]) -- post,pre,point
                end
                uprev = ucurrent
            end
            for i=1,#p[sp] do
                pcurrent = p[sp][i]
                if not pprev then
                    ppath = ppath .. string.format("M %d,%d",pcurrent[1],pcurrent[2])
                else
                    ppath = ppath .. string.format(" C %d,%d %d,%d %d,%d",
                        pprev[5],pprev[6],pcurrent[3],pcurrent[4],pcurrent[1],pcurrent[2]) -- post,pre,point
                end
                pprev = pcurrent
            end

            if sp == nofsp then close_path = " Z" else close_path = " Z " end
            upath,ppath = upath .. close_path, ppath .. close_path
        end


        return upath,ppath, bbx,bby, shiftx,shifty

    else
        return nil
    end
end


-- the path is saved as *.mvg file
local function save_IM_path(path,filename)
    if filename then
        -- check directory structure
        local absolute_shadowpath = lfs.currentdir() .. dir_separator .. shadowdir
        ensure_dir(absolute_shadowpath,shadowdir)
        local pathfile = absolute_shadowpath .. dir_separator .. filename

        local impath, err = io.open(pathfile,"w")
        if impath then
            impath:write("fill-rule nonzero path \'" .. path .. "\'")
            impath:close()
            job.variables.tobesaved.DropsPathChecksum =
            file.checksum(format("%s/path-%s.mvg",shadowdir,filename)) or "no checksum"
        else
            showmessage("drops","ioopenerror",{pathfile,err})
        end
    end
end



function drops.shadow(specification)
    local spec = specification

    local width,height,radius,mppath,resolution
    local penumbra,psigma,pdistance
    local    umbra,usigma,udistance

    local filecounter = drops.filecounter

    -- already in pixel (checked on TeX side)
    psigma = tonumber(spec.psigma)
    usigma = tonumber(spec.usigma)
    pdistance = tonumber(spec.pdistance) or fallback.pdistance
    udistance = tonumber(spec.udistance) or fallback.udistance
    resolution = tonumber(spec.resolution)

    -- the dimensions for the boxshadow graphic need to be accurate, to guarantee the same mppath name
    width = numberofpixels(spec.width)
    height = numberofpixels(spec.height)
    radius = numberofpixels(spec.radius) -- boxshadow only

    local path_id
    mppath = spec.mppath

    -- if no path is given, the default boxshadow template snaps in..
    if mppath == "" then
        mppath = drops.getboxshadowid(spec.width,spec.height,spec.radius)
        drops.createboxshadowpath(spec.width,spec.height,spec.radius)
        pathid = format("%s_%dppi",mppath,resolution)
        if not drops.IMpaths[pathid] then
            context(adapt_path_to_IM(mppath,pathid)) -- no scaling needed
        end
    else
        path_id = format("%s_%dppi",mppath,resolution)
    end

    -- check canvas sizes of the subshadows;
    -- the final IM canvas is always derived from the bigger penumbra shadow
    if udistance > pdistance then
        udistance = tonumber(fallback.udistance) ; pdistance = tonumber(fallback.pdistance)
        showmessage("drops","wrongcanvas")
    end

    -- percent (needs parameter check)
         umbra = tonumber(spec.umbra)    or tonumber(fallback.umbra)
      penumbra = tonumber(spec.penumbra) or tonumber(fallback.penumbra)

    if umbra <0 or umbra >100 then
        showmessage("drops","wrongumbra",umbra)
        umbra = tonumber(fallback.umbra)
    end
    if penumbra <0 or penumbra >100 then
        showmessage("drops","wrongpenumbra",penumbra)
        penumbra = tonumber(fallback.penumbra)
    end

    local fileformat,colorspace,shadowcolor,shadowbackgroundcolor
               fileformat = lower(spec.fileformat       or fallback.fileformat)
              shadowcolor = spec.shadowcolor            or fallback.shadowcolor
    shadowbackgroundcolor = spec.shadowbackgroundcolor  or fallback.shadowbackgroundcolor -- only relevant for flattening  (JPG)

    local color2text = drops.color2text
    -- a declared color space is the sign here, that the given color is handled as text; passed to IM as it is
    -- this way older ConTeXt versions can also use any color
    if spec.colorspace == "" then
        -- convert ConTeXt colors to text form, e.g. "red" -> "rgb(255,0,0)"; also derives the used colorspace from the given color
        shadowcolor, colorspace = color2text(shadowcolor)
          shadowbackgroundcolor = color2text(shadowbackgroundcolor)
    else
        colorspace = spec.colorspace -- color is given in text form (X11 color name or as 'rgb(),gray(),cmyk()');
    end

    if not (fileformat == "png" or fileformat == "jpg") then
        showmessage("drops","wrongfileformat",fileformat)
        fileformat = fallback.fileformat
    end
    if not (colorspace == "rgb" or colorspace == "cmyk" or colorspace == "gray") then
        showmessage("drops","wrongcolorspace",colorspace)
        colorspace = fallback.colorspace
    end
    if fileformat == "png" and colorspace == "cmyk" then
        showmessage("drops","pngandcmyk")
        fileformat = "jpg"
    end


    if trace_drops then
        report_drops("drops.boxshadow: w=%s, h=%s, r=%s, d=%s, umbra=%s, penumbra=%s, usigma=%s, psigma=%s, udistance=%d, pdistance=%d, col=%s, bgcol=%s, format=%s, colspace=%s, shbgcol=%s, mppath=%s",width,height,radius,resolution,umbra,penumbra,usigma,psigma,udistance,pdistance,shadowcolor,shadowbackgroundcolor,fileformat,colorspace,shadowbackgroundcolor,mppath)
    end


    -- internal shadow id
    local shadow_id = format("%s-%s-%s_%s_o%d:%d_s%d:%s_d%d:%s_%dppi_%s",
                       colorspace,shadowcolor,shadowbackgroundcolor,mppath,umbra,penumbra,usigma,psigma,udistance,pdistance,resolution,fileformat)
    current_shadow_id = shadow_id

    -- store the parameter for later usage (mainly for creating presets and also testing)
    local function store_current_parameters()
        local c = drops.parameters.current
        c.shadow_id = shadow_id; c.shadow_name = shadowname
        c.umbra = umbra; c.penumbra = penumbra
        c.usigma = usigma; c.psigma = psigma
        c.udistance = udistance; c.pdistance = pdistance
        c.resolution = resolution
        c.shadowcolor = shadowcolor; c.colorspace = colorspace
        c.shadowbackgroundcolor = shadowbackgroundcolor
    end


    local filename_template = "%s%d.%s" -- prefix number fileformat
    local shadowname = ""
    shadowsubdir = tex.jobname
    local absolute_shadowpath = lfs.currentdir() .. dir_separator .. shadowdir
    local absolute_shadowsubpath = absolute_shadowpath  .. dir_separator .. shadowsubdir

    if handled[shadow_id] then return end -- reusage of a new graphic, all done in a prior step
    -- shadow hash
    if drops.shadows[shadow_id] then
        shadowname = drops.shadows[shadow_id].name
        current_shadow_name = shadowname
        drops.parameters.current.shiftx = drops.shadows[shadow_id].shiftx or 0
        drops.parameters.current.shifty = drops.shadows[shadow_id].shifty or 0
        -- all parameter are set now, so we can store some data for later usage
        store_current_parameters()
        if trace_drops then report_drops("hash hit [" .. shadow_id .. "]= " .. shadowname) end

        if io.exists(absolute_shadowpath .. dir_separator .. im_batchfile) then -- are we ready for the check (batch was started)?
            if io.exists(absolute_shadowsubpath .. dir_separator .. shadowname) then -- check existence of gfx file
                return -- all ok
            else -- entry in hash, but no gfx file --> start recreation
                report_drops(format("warning: file '%s' for ID '%s' not found",shadowname,shadow_id))
                filecounter = 0 -- set number for batch comment
            end
        else return -- not ready yet
        end
    else
        filecounter = filecounter + 1
        drops.filecounter = filecounter
        shadowname = format(filename_template,
                            shadowfileprefix,filecounter,fileformat)
        drops.shadows[shadow_id] = { name = shadowname }
        current_shadow_name = shadowname

        -- all parameter are set now, so we can store some data for later usage
        store_current_parameters()
        if trace_drops then report_drops("new hash entry [" .. shadow_id .. "]= " .. shadowname) end
    end

    -- new graphic OR the graphic is already in the cache, but the file is not existent anymore
    run_batch_flag = true -- set the trigger to start batch file at end of current run
    handled[shadow_id] = true -- and mark graphic as handled


    -- IM batch file creation starts here...

    -- some IM related version handling
    local im_version = drops.im_version
    local IM_command, setcolorspace, cspace

    -- support manual setup with environment variable; maybe handy for local installations
    IM_command = os.getenv("IMCONV")
    if not IM_command or IM_command == "" then
        IM_command = im_version >= 7000000 and "magick" or "convert" end

    -- thanks to Fred Weinhaus
    -- RGB in combination with JPG is always sRGB (JPG does not support linear RGB)
    setcolorspace = ( im_version < 6070607 or im_version > 6070707 ) and " -set colorspace rgb" or ""
    cspace        = ( im_version < 6070606 or im_version > 6070707 ) and "rgb" or "srgb"
    if colorspace == "rgb" then colorspace = cspace end

    -- some abbreviations
    local lp,rp,cl,dq = left_parenthesis,right_parenthesis,continue_line,double_quote

    local quoted_shadowcolor = dq .. shadowcolor .. dq
    local quoted_shadowbackgroundcolor = dq .. shadowbackgroundcolor .. dq

    -- subdir for the graphics
    local relative_shadowpath, relative_shadowname
    relative_shadowpath = "."                 .. dir_separator .. shadowdir .. dir_separator .. shadowsubdir
    relative_shadowname = relative_shadowpath .. dir_separator .. shadowname

    -- set the density for the saved graphic by using both 'units' and 'density'
    -- '-set density' is ignored by IM 6.6.0-4 (Debian)
    local density_options = " -units PixelsPerInch -density " .. resolution

    -- set the final color space
    local colorspace_options = " -colorspace " .. colorspace .. " "
    -- PNG
    -- '-colorspace' must be last option (otherwise problems with CMYK and shadowcolor=black)!
    local png_color_options = " -background " .. quoted_shadowbackgroundcolor .. colorspace_options
    -- first digit is zlib compression level, second stands for encoding filter (see http://www.imagemagick.org/Usage/formats/#png for details)
    -- no need to waste time here, as luatex will recompress the graphic in any case;
    local png_quality_options = " -quality 00"
    -- no IM date stamps ('tEXt' chunk: 2*37bytes for creation|modify date)
    -- the png 'tIME' chunk (7bytes), which contains the time of the last modification, is still written
    local png_no_IM_time_stamps = " +set date:create +set date:modify"
    -- JPG
    -- '-set colorspace' is needed to guarantee a linear working color space for the conversion (else we get wrong colors)
    -- IM 6.7.7-8 and newer: " -alpha Remove"
    -- IM 6.6.0-4: ' -flatten'
    local jpg_color_options = setcolorspace .. " -background " .. quoted_shadowbackgroundcolor .. " -alpha Remove" .. colorspace_options
    local jpg_quality = 95 -- noticeable compression artefacts even in the lower nineties...
    local jpg_quality_options = " -quality " .. jpg_quality


    local gfx_options
    if fileformat == "png" then
        -- '-colorspace' must be last parameter (otherwise problems with CMYK and shadowcolor=black)!
        gfx_options = density_options .. png_quality_options .. png_no_IM_time_stamps  .. png_color_options
    else
        gfx_options = density_options .. jpg_quality_options .. jpg_color_options
    end

    -- only needed for the fallback code (path creation fails)
    local border = 2 -- add x pixels around mask to get reliable corners (can differ with small radii otherwise)

    -- '-depth 8' is needed for the correct interpretation of the shadow color values
    -- must stand before 'xc:none'; using the '--verbose' parameter showed 16-bit for the first shadow mask otherwise
    local mask_setup = [[-size %dx%d -depth 8 xc:none -fill black -stroke none -draw ]]
--    local mask_setup = [[-size %dx%d -depth 8 xc:none -fill none -stroke white -draw ]] -- debug
    local mask_fallback = [["roundrectangle 1,1 %d,%d %d,%d"]]

    local mask_template
    -- shifting the shadow in the graphic makes no sense here, as the positioning is done by ConTeXt;
    local shadow_template = [[-background %s -shadow %dx%d+0+0 +repage]]


    batch[#batch+1] = format("\n%s [%d]  %s",
                                 comment_sign,filecounter,mppath)

    -- calculate the paths for umbra and penumbra masks
    local umask_path,pmask_path,umbra_mask,penumbra_mask
    local twidth,theight
    umask_path,pmask_path,twidth,theight,shiftx,shifty = get_IM_paths(mppath,udistance,pdistance)

    if (umask_path and pmask_path) then
        width,height = twidth,theight
    else
        showmessage("drops","noshadowpath")
    end

    --  1.
    if pmask_path then
        if #pmask_path > IM_MAX_PATH_LIMIT then
            local mask_path_filename = format("path-%s-%s-%sppi-pd%i.mvg",
                                              tex.jobname,mppath,resolution,pdistance)
            save_IM_path(pmask_path,mask_path_filename)
            mask_template = mask_setup .. "@" .. shadowdir .. dir_separator ..  mask_path_filename
        else
            mask_template = mask_setup .. "\"fill-rule nonzero path \'" .. pmask_path .. "\' \""
        end
        -- store the shift values of the bigger penumbra shadow
        drops.shadows[shadow_id].shiftx = shiftx ; drops.shadows[shadow_id].shifty = shifty
        drops.parameters.current.shiftx = shiftx ; drops.parameters.current.shifty = shifty
        penumbra_mask = format(mask_template,
                                      width,height)
    else -- fallback in case of a path problem
        mask_template = mask_setup .. mask_fallback
        penumbra_mask = format(mask_template,
                                      width+2*pdistance+border,height+2*pdistance+border, -- total size
                                      width+2*pdistance,height+2*pdistance, radius,radius) -- drawn mask size
    end

    -- 2.
    local penumbra_shadow = format(shadow_template,
                                   quoted_shadowcolor,penumbra,psigma)
    -- 3.
    if umask_path then
        if #umask_path > IM_MAX_PATH_LIMIT then
            local mask_path_filename = format("path-%s-%s-%sppi-ud%i.mvg",
                                              tex.jobname,mppath,resolution,udistance)
            save_IM_path(umask_path,mask_path_filename)
            mask_template = mask_setup .. "@" .. shadowdir .. dir_separator ..  mask_path_filename
        else
            mask_template = mask_setup .. "\"fill-rule nonzero path \'" .. umask_path .. "\' \""
        end
        umbra_mask = format(mask_template,
                                   width,height)

    else -- fallback in case of a path problem
        mask_template = mask_setup .. mask_fallback
        umbra_mask = format(mask_template,
                                   width+2*udistance+border,height+2*udistance+border, -- total size
                                   width+2*udistance,height+2*udistance, radius,radius) -- drawn mask size
    end

    -- .4
    local umbra_shadow = format(shadow_template,
                                quoted_shadowcolor,umbra,usigma)
    -- 5.
    -- does not work in 6.6.0-4 (Debian)
    -- IM 6.9.2-0 needs background color
    local combine_shadows =
        format(" -background %s -channel Alpha -gravity Center -compose Lighten -composite",quoted_shadowcolor)
    -- 6
    local blur_radius,blur_sigma = 0, 2
    local frame_size = blur_sigma+1
    local optimize_shadow =
--        format("+clone -bordercolor none -border %d -background %s -alpha Background -channel Alpha -blur %dx%.1f", -- IM 6.9.2-0: problem with '-border'
        format("+clone -bordercolor none -compose Src -frame %dx%d -background %s -alpha Background -channel Alpha -blur %dx%.1f",
               frame_size,frame_size,quoted_shadowcolor,blur_radius,blur_sigma)

    -- IM SHADOW WITH UMBRA + PENUMBRA AREA
    --  1. create penumbra mask as base for shadow
    --  2. create penumbra shadow
    --  3. create umbra mask
    --  4. create umbra shadow
    --  5. combine transparent channels of both shadows
    --  6. adds full transparent border (simply avoids PDF render problems at graphic borders),
    --     corrects the border color (black->shadowcolor) afterwards to save some bytes
    --     and also blurs possibly visible transitions between the sub-shadow areas
    --  7. save shadow

    batch[#batch+1] = IM_command .. cl ..
        lp .. penumbra_mask                .. rp .. cl .. -- (1.)
        lp .. "+clone " .. penumbra_shadow .. rp .. cl .. -- (2.)
        lp .. umbra_mask                   .. rp .. cl .. -- (3.)
        lp .. "+clone " .. umbra_shadow    .. rp .. "-delete 0,2 " .. cl .. -- (4.)
        combine_shadows                    .. cl .. -- (5.)
        lp .. optimize_shadow              .. rp .. "-delete 0 "   .. cl .. -- (6.)
        gfx_options .. dq .. relative_shadowname .. dq -- (7.)



end
