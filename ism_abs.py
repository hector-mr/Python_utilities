def mm_abs(coll_den, Energy = []):
    """
    Calculate column density absorption in X-ray wavelengths.
    The reference is Morrison & McCammon 1983 ApJ, 270, 119
    http://adsabs.harvard.edu/cgi-bin/nph-basic_connect
    """
 
    #------------------------------------------------------------------------------------------------------------------------------#
    
    # When you define your function using this syntax:

    # def someFunc(*args)
    # for x in args
    #     print x

    # You are telling it that you expect a variable number of arguments. 
    # If you want to pass in a List (Array from other languages) you'd do something like this:

    # def someFunc(myList = [], *args)
    #     for x in myList:
    #         print x

    # Then you can call it with this:

    # items = [1,2A,3,4,5]

    # someFunc(items)

    # You need to define named arguments before variable arguments, and variable arguments before keyword arguments. 
    # You can also have this:

    # def someFunc(arg1, arg2, arg3, *args, **kwargs)
    #     for x in args
    #         print x


    #------------------------------------------------------------------------------------------------------------------------------#
    
    
    import numpy as np
    
    # Beware of numpy.float32/64
    if (type(Energy) is float or type(Energy) is np.float32 or type(Energy) is np.float64 or type(Energy) is int) : Energy = [Energy]
        
        
    sigma = []
    photar = []
        

        
    for idx, E in enumerate(Energy):
        
        
        if ((E < 0.03) | (E > 10.0)) :     
            c0 = 0.0
            c1 = 0.0
            c2 = 0.0
    
        elif ((E >= 0.03) & (E <= 0.100)) :
            c0 = 17.3
            c1 = 608.1
            c2 = -2150.0
        
        
        elif ((E > 0.100) & (E <= 0.284)) :
            c0 = 34.6
            c1 = 267.9
            c2 = -476.1
        
        elif ((E > 0.284) & (E <= 0.400)) :
            c0 = 78.1
            c1 = 18.8
            c2 = 4.3
        
        elif ((E > 0.400) & (E <= 0.532)) :
            c0 = 71.4
            c1 = 66.8
            c2 = -51.4
        
        elif ((E > 0.532) & (E <= 0.707)) :
            c0 = 95.5
            c1 = 145.8
            c2 = -61.1
        
        elif ((E > 0.707) & (E <= 0.867)) :
            c0 = 308.9
            c1 = -380.6
            c2 = 294.0
        
        elif ((E > 0.867) & (E <= 1.303)) :
            c0 = 120.6
            c1 = 169.3
            c2 = -47.7
        
        elif ((E > 1.303) & (E <= 1.840)) :
            c0 = 141.3
            c1 = 146.8
            c2 = -31.5
        
        elif ((E > 1.840) & (E <= 2.471)) :
            c0 = 202.7
            c1 = 104.7
            c2 = -17.0
        
        elif ((E > 2.471) & (E <= 3.210)) :
            c0 = 342.7
            c1 = 18.7
            c2 = 0.0
        
        elif ((E > 3.210) & (E <= 4.038)) :
            c0 = 352.2
            c1 = 18.7
            c2 = 0.0
        
        elif ((E > 4.038) & (E <= 7.111)) :
            c0 = 433.9
            c1 = -2.4
            c2 = 0.75
        
        elif ((E > 7.111) & (E <= 8.331)) :
            c0 = 629.0
            c1 = 30.9
            c2 = 0.00
        
        elif ((E > 8.331) & (E <= 10.00)) :
            c0 = 701.2
            c1 = 25.2
            c2 = 0.00
        
        else:
            c0 = 0.0
            c1 = 0.0
            c2 = 0.0
        
        
        sigma.append((c0 + c1*E + c2*E**2)/E**3)

        # "coll_den" is the nH column density in cm2

        # sigma is in units of 1E-24 cm2
        # the absorption is given as A(E) = EXP(-coll_den*sigma(E))
    
    
        # For really small numbers
        photar.append(np.maximum(1.E-10, np.exp(-coll_den*sigma[idx]*1.0E-24)))
    
    
    photar = np.asarray(photar)
    return photar