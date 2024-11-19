def GetStrawCenters(strawRadius, verticalOffset): # All in cm units
    # Define the straws to cover +/- 10cm horizontally, in three staggered layers
    # Plot to confirm the geometry too
    yCenter, zCenter = [], []
    nTubes = int(5/strawRadius)+1;
    # Generate the first layer
    for i in range(2*nTubes+1):
        yCenter.append(verticalOffset)
        zCenter.append(-nTubes*2*strawRadius + i*2*strawRadius)
    # Generate the second layer
    for i in range(2*nTubes+1):
        yCenter.append(verticalOffset+np.sqrt(3)*strawRadius)
        zCenter.append(-nTubes*2*strawRadius + i*2*strawRadius + strawRadius)
     # Generate the third layer
    for i in range(2*nTubes+1):
        yCenter.append(verticalOffset + 2*np.sqrt(3)*strawRadius)
        zCenter.append(-nTubes*2*strawRadius + i*2*strawRadius)
    return zCenter, yCenter

def PlotStraws(zCenter, yCenter, radius):
    figure, axes = plt.subplots()
    axes.set_aspect(1)
    axes.set_xlim([-15,15])
    axes.set_ylim([0,30])
    for i in range(len(yCenter)):
        circleTube = plt.Circle((zCenter[i], yCenter[i]), radius, fill = False)
        axes.add_artist(circleTube)
    plt.show()

def PlotStrawsAndTrack(zCenter, yCenter, radius, hitI, hitR, m, zOrigin):
    figure, axes = plt.subplots()
    axes.set_aspect(1)
    axes.set_xlim([-15,15])
    axes.set_ylim([0,30])
    for i in range(len(yCenter)):
        circleTube = plt.Circle((zCenter[i], yCenter[i]), radius, fill = False)
        axes.add_artist(circleTube)
        for j in range(len(hitI)):
            if i==hitI[j]:
                trackCircle = plt.Circle((zCenter[i], yCenter[i]), hitR[j], fill = False, color="red")
                axes.add_artist(trackCircle)
    p1=[zOrigin,0]
    lastIndex = len(yCenter)-1
    p2=[(yCenter[lastIndex]+m*zOrigin)/m,yCenter[lastIndex]]
    plt.axline(p1, p2, color="red")
    title = "Plot and Trajectory for rTube:{:.2f}".format(strawRadius)
    plt.title(title)
    filename = "StrawsandTrack_rTube_{}.png".format(file_number)
    plt.savefig(save_results_to + filename, format= "png")
    plt.show()


def GetRadiiForTrack(zCenter, yCenter, radius, m, zOrigin):
    b = -m*zOrigin
    # Loop over all the points and find the tubes we went through (distance to the tube center < radius)
    hitI, hitR = [], []
    nHit = 0
    for i in range(len(yCenter)):
        d = np.abs(m*zCenter[i] - yCenter[i] + b)/(np.sqrt(m*m + 1))
        if (d <= radius):
            hitI.append(i)
            hitR.append(d)
            nHit+=1
    #print("Tubes hit = ", nHit)
    return nHit, hitI, hitR



def fLinear(x,A,B):
    return x*A+B


def FitTrackFromRadii(zCenter, yCenter, hitI, hitR, strawR):
    # The approach here is to minimize the difference between the 'measured' r and the fit track r for all
    # the hit tubes to constrain the track direction/source
    # Figure out a guess for slope from center of the hit tubes
    xVal, xValNew, yVal = [], [], []
    for i in range(len(hitI)):
        xVal.append(zCenter[hitI[i]])
        xValNew.append(0)
        yVal.append(yCenter[hitI[i]])
    popt, pcov = opt.curve_fit(fLinear, xVal, yVal)
    m0 = popt[0]
    b = popt[1]
    zOff0 = -b/m0
    if (np.abs(zOff0) > 10):
        zOff0 = 0
    fomToBeat = CalcFOM([m0, zOff0], zCenter, yCenter, hitI, hitR)
    # Now iterate through fit lines going through at +/- the hitR of the straw tube and find the closest
    for i0 in range(-1, 2, 2):
        for i1 in range(-1, 2, 2):
            for i2 in range(-1, 2, 2):
                for i3 in range(-1, 2, 2):
                    for i4 in range(-1, 2, 2):
                        xValNew[0] = xVal[0] + i0*hitR[0]
                        if (len(hitI) > 1):
                            xValNew[1] = xVal[1] + i1*hitR[1]
                        if (len(hitI) > 2):
                            xValNew[2] = xVal[2] + i2*hitR[2]
                        if (len(hitI) > 3):
                            xValNew[3] = xVal[3] + i3*hitR[3]
                        if (len(hitI) > 4):
                            xValNew[4] = xVal[4] + i4*hitR[4]
                        popt, pcov = opt.curve_fit(fLinear, xValNew, yVal)
                        m = popt[0]
                        b = popt[1]
                        zOff = -b/m
                        fom = CalcFOM([m, zOff], zCenter, yCenter, hitI, hitR)
                        if (fom < fomToBeat):
                            m0 = m
                            zOff0 = zOff
                            fomToBeat = fom

    # And finetune now, though the fomToBeat should be pretty good
    result = opt.minimize(CalcFOM,
                            x0=[m0,zOff0],
                            args=(zCenter, yCenter, hitI, hitR),
                            method='Nelder-Mead')

    return result.x[0], result.x[1]


def CalcFOM(x, zCenter, yCenter, hitI, hitR):
    # do lines match with observed radii
    m = x[0]
    zOrigin = x[1]
    fom = 0
    b = -m*zOrigin
    for i in range(len(hitI)):
        d = np.abs(m*zCenter[hitI[i]] - yCenter[hitI[i]] + b)/(np.sqrt(m*m + 1))
        fom += np.abs(hitR[i]-d)
    #if (zOrigin < -10 or zOrigin > 10):
     #   fom += 1000.
    return fom



#define gaussian for histogram centerpoint scatter curvefit
def Gauss(x, a, avg, sigma):
    return a*np.exp(-(x-avg)**2/(2*sigma**2))



#use to do range in for-loops with float steps rather than integer iterations
def range_with_floats(start, stop, step):
    while stop > start:
        yield start
        start += step

def CalculateSignalWidth(strawRadius, pathRadius):
    # Here we want to insert the GARFIELD code to calculate a signal at the given radius and extract the signal width

    # Getting our specific gas
    gas = ROOT.Garfield.MediumMagboltz()
    gas.LoadGasFile('ar_93_co2_7_3bar.gas')
    gas.LoadIonMobility(path + '/share/Garfield/Data/IonMobility_Ar+_Ar.txt')


    cmp = ROOT.Garfield.ComponentAnalyticField()
    cmp.SetMedium(gas)

    # Wire radius [cm]
    rWire = 25.e-4



    # Outer radius of the tube [cm]
    #strawRadius = 0.75



    # Voltages (min:1000, max:3000)
    vWire = 2730.
    vTube = 0.
    # Add the wire in the centre.
    cmp.AddWire(0., 0., 2 * rWire, vWire, 's')
    # Add the tube.
    cmp.AddTube(strawRadius, vTube, 0)


    # Make a sensor.
    sensor = ROOT.Garfield.Sensor()
    sensor.AddComponent(cmp);
    sensor.AddElectrode(cmp, 's')
    # Set the signal time window.
    tstep = 0.5;
    tmin = -0.5 * tstep
    nbins = 1000
    sensor.SetTimeWindow(tmin, tstep, nbins)
    # Set the delta reponse function.
    infile = open('mdt_elx_delta.txt', 'r')
    times = ROOT.std.vector('double')()
    values = ROOT.std.vector('double')()
    for line in infile:
      line = line.strip()
      line = line.split()
      times.push_back(1.e3 * float(line[0]))
      values.push_back(float(line[1]))
    infile.close()
    sensor.SetTransferFunction(times, values)



    # Set up Heed.
    track = ROOT.Garfield.TrackHeed()
    track.SetParticle('proton')
    #around 100 MeV
    track.SetKineticEnergy(100e6)
    track.SetSensor(sensor)



    # RKF integration.
    drift = ROOT.Garfield.DriftLineRKF()
    drift.SetSensor(sensor)
    drift.SetGainFluctuationsPolya(0., 20000.)
    # drift.EnableIonTail()


    driftView = ROOT.Garfield.ViewDrift()
    cD = ROOT.TCanvas('cD', '', 600, 600)
    driftView.SetCanvas(cD)
    cellView = ROOT.Garfield.ViewCell()
    plotDrift = False
    if plotDrift:
      drift.EnablePlotting(driftView)
      track.EnablePlotting(driftView)
      cellView.SetComponent(cmp)
      cellView.SetCanvas(driftView.GetCanvas())



    plotSignals = True



    x0 = pathRadius
    y0 = -math.sqrt(strawRadius * strawRadius - pathRadius * pathRadius)

    y_basket = []
    x_basket = []
    width = 0



    nTracks = 1
    for j in range(nTracks):
      sensor.ClearSignal()
      track.NewTrack(x0, y0, 0, 0, 0, 1, 0)
      for cluster in track.GetClusters():
        for electron in cluster.electrons:
          drift.DriftElectron(electron.x, electron.y, electron.z, electron.t)
        if plotDrift:
          driftView.GetCanvas().Clear()
          cellView.Plot2d()
          driftView.Plot(True, False)
          ROOT.gPad.Update()
      sensor.ConvoluteSignals()
      nt = ctypes.c_int(0)
      if sensor.ComputeThresholdCrossings(-2., 's', nt) == False:
        continue
      if plotSignals:
        # Matplotlib plotting (Python stuff)
        x = np.arange(-0.5*tstep, (nbins-1)*tstep, tstep)
        y = np.zeros(1000)
        for i in range(1000):
            y[i] = sensor.GetSignal('s',i)

            if y[i] <= -0.25:
                #y_basket.append(y[i])
                x_basket.append(x[i])
                continue

        if len(x_basket) >= 2:
            width = x_basket[-1] - x_basket[0]

        else:
            width = 0

    signalWidth = width


    return signalWidth

def GetRadiusFromWidth(strawRadius,signalWidth):
    # Here we want to use the relationship extracted for average signal width to radius for a given strawtube size


    if strawRadius == 0.50:
        slope = -0.005996070857
        intercept = 0.6848588078

    if strawRadius == 0.55:
        slope = -0.005756078266
        intercept = 0.7269435018

    if strawRadius == 0.60:
        slope = -0.005580029419
        intercept = 0.7765892342

    if strawRadius == 0.65:
        slope = -0.005118370892
        intercept = 0.8159194578

    if strawRadius == 0.70:
        slope = -0.0052958385
        intercept = 0.8858135288

    if strawRadius == 0.75:
        slope = -0.005074717041
        intercept = 0.9303558177

    if strawRadius == 0.80:
        slope = -0.004914084831
        intercept = 0.9763605203

    if strawRadius == 0.85:
        slope = -0.004579289866
        intercept = 1.019672194

    if strawRadius == 0.90:
        slope = -0.004300621838
        intercept = 1.053465811

    if strawRadius == 0.95:
        slope = -0.00416347126
        intercept = 1.104273789

    if strawRadius == 1.00:
        slope = -0.003808395786
        intercept = 1.139277676

    if strawRadius == 1.05:
        slope = -0.003470040325
        intercept = 1.174006702

    if strawRadius == 1.10:
        slope = -0.003180942378
        intercept = 1.210846744

    if strawRadius == 1.15:
        slope = -0.002675870475
        intercept = 1.227340657



    rTrack = slope * signalWidth + intercept
    #print(rTrack)

    return rTrack

  eighty_percent_list = [3.5000000000000018, 6.699999999999992, 5.999999999999995, 7.19999999999999, 7.39999999999999, 4.799999999999999, 4.899999999999999, 4.799999999999999, 5.4999999999999964, 6.8999999999999915, 4.5, 5.999999999999995, 5.699999999999996]

rTube_list = [1.15]
event_number = 500


for rTube in rTube_list:



    # Setup our straw tube layers

    strawRadius = rTube   # 2 places past decimal, ranges from 0.5 - 1.15, in .05 increments (in cm)


    # Filename Conversion
    #number = str(strawRadius*100)
    file_number = "{:0>{}}".format(strawRadius, 3)




    zStraw, yStraw = GetStrawCenters(strawRadius, 10.)

    # Calculate a bunch of trajectories for a proton from an initial position at z=zPos, y=0, going in
    # a random direction toward the tubes.  For each, we calculate the true radius in each tube hit, then
    # calculate a signal to get a width, and relate this back to the radius based on the larger set of simulations.
    # We then recalculate the track and see the spread we end up with in determining zPos.

    zPos = 0.
    zVal = []
    pathRadiuslist = []
    newTrackRadiilist = []

    for trackNum in range(event_number):
        nHit = 0
        while (nHit < 2):
            slope = rndm.uniform(-10,10)
            nHit, index, pathRadius = GetRadiiForTrack(zStraw, yStraw, strawRadius, slope, zPos)

        newTrackRadii = []

        for i in range(nHit):
            tWidth = CalculateSignalWidth(strawRadius, pathRadius[i])
            pathRadiuslist.append(pathRadius[i])
            #print(GetRadiusFromWidth(strawRadius, tWidth))
            newTrackRadii.append(GetRadiusFromWidth(strawRadius, tWidth))
            newTrackRadiilist.append(newTrackRadii[i])


        m, zFit = FitTrackFromRadii(zStraw, yStraw, index, newTrackRadii, strawRadius)
        zVal.append(zFit)

    # CHECK VERSION
    #     m, zFit = FitTrackFromRadii(zStraw, yStraw, index, pathRadius, strawRadius)
    #     zVal.append(zFit)




    # Histogram of counts
    n, bins, patches = plt.hist(x=zVal, bins=200, range=(zPos-20, zPos+20), color='#0504aa',
                                alpha=0.7, rwidth=0.85)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('Z Value')
    plt.ylabel('Counts')
    maxfreq = n.max()
    # Set a clean upper y-axis limit.
    plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
    plt.xlim(zPos-20, zPos+20)
    title = "Histogram of Z-values Counts for rTube:{:.2f}".format(strawRadius)
    plt.title(title)
    filename = "Histogram_rTube_{}.png".format(file_number)
    plt.savefig(save_results_to + filename, format= "png")
    plt.show()





    # Trajectory Diagram
    nHit, index, pathRadius = GetRadiiForTrack(zStraw, yStraw, strawRadius, slope, zPos)
    m, zFit = FitTrackFromRadii(zStraw, yStraw, index, pathRadius, strawRadius)
    PlotStrawsAndTrack(zStraw, yStraw, strawRadius, index, pathRadius, slope, zPos)






    #Gaussian CurveFit of center points of the initial histogram bins
    x = bins[:-1]
    y = n

    # Plot out the current state of the data and model
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x, y)

    # popt are parameters of Gauss function, pcov is the covariance
    popt, pcov = curve_fit(Gauss, x, y)

    #popt returns the best fit values for parameters of the given model (func)
    print(' "a" coefficient =', popt[0]," mean =", popt[1]," std =", popt[2])

    y_fit = Gauss(x, popt[0], popt[1], popt[2])
    ax.plot(x, y_fit, c='r', label='Best fit')
    plt.xlabel('Z Value')
    plt.ylabel('Counts')
    maxfreq = n.max()
    # Set a clean upper y-axis limit.
    plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
    plt.xlim(zPos-20, zPos+20)
    ax.legend()
    title = "Gaussian Centerpoint Curvefit for rTube:{:.2f}".format(strawRadius)
    plt.title(title)
    filename = "GaussianCurvefit_rTube_{}.png".format(file_number)
    plt.savefig(save_results_to + filename, format= "png")
    plt.show()







    # Calculated Info from Scatter Plot and Curve Fit
    # Given information
    avg = popt[1]
    std = popt[2]

    # Calculate z-score for min
    z_min = avg - (3 * std)
    # Calculate z-score for max
    z_max = avg + (3 * std)

    print("Lower End", z_min)
    print("Upper End", z_max)

    n_in_range = 0

    for i in zVal:
        if i < z_max and i > z_min:
            n_in_range += 1

    #print("n_in_range=", n_in_range)

    percent = (n_in_range / len(zVal)) * 100
    print("Percentage of counts within 3 standard devaitions is", round(percent, 2), "%")


    for i in range_with_floats(0, 10, 0.1):

        n_range = 0
        #print(i)
        data_max = avg + i
        data_min = avg - i
        for j in zVal:
            if j < data_max and j > data_min:
                n_range += 1
        #print(n_range)

        if (n_range / len(zVal)) > 0.8:
            width80 = i
            #print(width80)
            #print("stop, this is the width that covers 80% of data")
            break

    std_amounts = (width80/std)/2
    print("80% of the data is covered within",u"\u00B1", width80 ,"of the of the zVal,\nor roughly", round(std_amounts, 2), "standard deviations out from average")
    eighty_percent_list.append(width80)









    # Finding difference between pathRadii and NewTrackRadii for each trial
    diff = []
    for i in range(len(pathRadiuslist)):
         diff.append(abs(pathRadiuslist[i] - newTrackRadiilist[i]))



    # CHECK: THESE VALUES SHOULD BE SAME
    # print(len(pathRadiuslist))
    # print(len(newTrackRadiilist))
    # print(len(diff))



    # Plotting absolute values of difference
    plt.figure()
    plt.xlim((0,1))
    plt.xlabel("Difference")
    plt.ylabel("Count")
    title = "Absolute Value Differece Between Path Radii (True) and newTrackRadii (Simulated) for rTube:{:.2f}".format(strawRadius)
    plt.title(title)
    plt.hist(diff,40)
    filename = "AbsoluteDiff_rTube_{}.png".format(file_number)
    plt.savefig(save_results_to + filename, format= "png")
    plt.show()




    # Plot absolute difference as a function of the pathRadius (scatterplot) [difference vs radius]
    plt.figure()
    #plt.xlim((0,len()))
    title = "Absolute Difference as a function of pathRadius for rTube:{:.2f}".format(strawRadius)
    plt.title(title)
    plt.xlabel("PathRadius")
    plt.ylabel("Absolute Difference from True Value")
    plt.plot(pathRadiuslist, diff,'bo')
    filename = "AbsoluteDiff_vs_PathRadius_rTube_{}.png".format(file_number)
    plt.savefig(save_results_to + filename, format= "png")
    plt.show()


print("eighty percent list so far:", eighty_percent_list)
