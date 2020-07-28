import numpy
import moments
from moments import Numerics, Manips, Integration
from moments.Spectrum_mod import Spectrum


######################
####    SI  and SI+   ####
######################

def SI(params, ns):
    """
    nu1= pop size for North Sea
	nu2=pop size for Baltic Sea 
	T1= time of split
    """
    nu1,nu2,T1 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01,m = numpy.array([[0, 0], [0, 0]]))
    return fs
    

   
def SI_b(params, ns):
    """
    nu1= pop size for North Sea
	s=proportion of the North Sea pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Baltic Sea pop
	T1= time of split
    """
    nu1,nu2,s,T1= params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate(nu_func, T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    return fs 


def SI_ae(params, ns):
    """
    nu1= pop size following ancestral population expnasion (assumed constant in both NS and Baltic Sea)
	Tae= timing of ancestral population expansion
	T1= time of split
    """
    nu_ae,nu1,nu2,Tae,T1= params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    return fs
    
    
def SI_ae_b(params, ns):
    """
    nu1= pop size after ancestral expansion (this remains constant for teh North sea population after the split)
	s=proportion of the North Sea pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Baltic Sea pop
	Tae= timing of ancestral population expansion
	T1= time of split and start of population growth in the Baltic Sea
    """
    nu_ae,nu1,nu2,s,Tae,T1= params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate(nu_func, T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    return fs


def SI_2N(params, ns):
    """
    nu1= pop size for North Sea
	nu2=pop size for Baltic Sea 
	T1= time of split
	hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
	Q= proportion of the genome under linked selection 
    """
    nu1,nu2,T1,hrf,Q = params
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
# calculate the spectrum for lower Ne regions of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2 = moments.Spectrum(stsn2)
    fsn2 = moments.Manips.split_1D_to_2D(fsn2, ns[0], ns[1])
    fsn2.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fs2= Q*fsn2+(1-Q)*fs
    return fs2
 
 
def SI_b_2N(params, ns):
    """
    nu1= pop size for North Sea
	s=proportion of the North Sea pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Baltic Sea pop
	T1= time of split
	hrf= Hill-Robertson factor (i.e. average Ne for regions under selection as a proportion of "neutral" Ne)
	Q= proportion of the genome under selection
    """
    nu1,nu2,s,T1,hrf,Q = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func= lambda t: [nu1,nu2_func(t)]
    nu2hrf_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)*hrf
    nuhrf_func= lambda t: [nu1*hrf,nu2hrf_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate(nu_func, T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
# calculate the spectrum from low Ne parts of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2 = moments.Spectrum(stsn2)
    fsn2 = moments.Manips.split_1D_to_2D(fsn2, ns[0], ns[1])
    fsn2.integrate(nuhrf_func, T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))

	
    fs2=Q*fsn2+(1-Q)*fs
    return fs2
    

def SI_ae_2N(params, ns):
    """
    nu1= pop size for North Sea
	Tae= time of ancestral pop expansion 
	T1= time of split
	hrf= Hill-Robertson factor (i.e. average Ne for regions under linked selection as a proportion of "neutral" Ne)
	Q= proportion of the genome under linked selection 
    """
    nu_ae,nu1,nu2,Tae,T1,hrf,Q = params
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
# calculate the spectrum for lower Ne regions of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2 = moments.Spectrum(stsn2)
    fsn2.integrate([nu1*hrf], Tae)
    fsn2 = moments.Manips.split_1D_to_2D(fsn2, ns[0], ns[1])
    fsn2.integrate([nu1*hrf, nu1*hrf], T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fs2= Q*fsn2+(1-Q)*fs
    return fs2




def SI_ae_b_2N(params, ns):
    """
    nu1= pop size for North Sea
	s=proportion of the North Sea pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Baltic Sea pop
	Tae= timing of ancestral population expansion
	T1= time of split
	hrf= Hill-Robertson factor (i.e. average Ne for regions under selection as a proportion of "neutral" Ne)
	Q= proportion of the genome under selection
    """
    nu_ae,nu1,nu2,s,Tae,T1,hrf,Q = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func= lambda t: [nu1,nu2_func(t)]
    nu2hrf_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)*hrf
    nuhrf_func= lambda t: [nu1*hrf,nu2hrf_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate(nu_func, T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
# calculate the spectrum from low Ne parts of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2 = moments.Spectrum(stsn2)
    fsn2.integrate([nu1*hrf], Tae)
    fsn2 = moments.Manips.split_1D_to_2D(fsn2, ns[0], ns[1])
    fsn2.integrate(nuhrf_func, T1, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))

	
    fs2=Q*fsn2+(1-Q)*fs
    return fs2

######################
####    IM  and IM+   ####
######################
    
def IM(params, ns):
    """
    nu1= pop size for North Sea
	nu2=pop size for Baltic Sea 
	T1= time of split
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
    """
    nu1,nu2,T1,m12,m21 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    return fs

def IM_2M(params, ns):
    """
	nu1= pop size for North Sea 
	nu2= pop size for the Baltic Sea
	= time of population split
	T1= time of split
	i1= reduction in migration rate for Islands from the Baltic to North Sea 
	i2= reduction in migration rate for Islands from the North Sea to Baltic
	P= proportion of the genome made up of "islands"
    """
    nu1,nu2,T1,m12,m21,i1,i2,P = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
# calculate the spectrum from genomic islands
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsi = moments.Spectrum(stsi)
    fsi = moments.Manips.split_1D_to_2D(fsi, ns[0], ns[1])
    fsi.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fs2=P*fsi+(1-P)*fs
    return fs2

def IM_ae(params, ns):
    """
    nu1= pop size following ancestral population expnasion (assumed constant in both NS and Baltic Sea)
	Tae= timing of ancestral population expansion
	T1= time of split
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
    """
    nu_ae,nu1,nu2,Tae,T1,m12,m21 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    return fs
def IM_ae_2M(params, ns):
    """
	nu1= pop size after ancestral expansion
	Tae= time of ancestral expansion
	T1= time of split
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
	i1= reduction in migration rate for Islands from the Baltic to North Sea 
	i2= reduction in migration rate for Islands from the North Sea to Baltic
	P= proportion of the genome made up of "islands"
    """
    nu_ae,nu1,nu2,Tae,T1,m12,m21,i1,i2,P = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
# calculate the spectrum from genomic islands
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsi = moments.Spectrum(stsi)
    fsi.integrate([nu_ae], Tae)
    fsi = moments.Manips.split_1D_to_2D(fsi, ns[0], ns[1])
    fsi.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fs2=P*fsi+(1-P)*fs
    return fs2

def IM_b(params, ns):
    """
    nu1= pop size for North Sea
	s=proportion of the North Sea pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Baltic Sea pop
	T1= time of split
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
    """
    nu1,nu2,s,T1,m12,m21 = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate(nu_func, T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    return fs

def IM_b_2M(params, ns):
    """
    nu1= pop size for North Sea
	s=proportion of the North Sea pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Baltic Sea pop
	= time of population split
	T1= time of second epoch of migration and population growth
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
	i1= reduction in migration rate for Islands from the Baltic to North Sea 
	i2= reduction in migration rate for Islands from the North Sea to Baltic
	P= proportion of the genome made up of "islands"
    """
    nu1,nu2,s,T1,m12,m21,i1,i2,P = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate(nu_func, T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
# calculate teh spectrum for genomic islands	
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsi = moments.Spectrum(stsi)
    fsi = moments.Manips.split_1D_to_2D(fsi, ns[0], ns[1])
    fsi.integrate(nu_func, T1, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
	
    fs2=P*fsi+(1-P)*fs
    return fs2



def IM_ae_b(params, ns):
    """
    nu1= pop size after ancestral expansion (this remains constant for teh North sea population after the split)
	s=proportion of the North Sea pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Baltic Sea pop
	Tae= timing of ancestral population expansion
	T1= time of split and start of population growth in the Baltic Sea
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
    """
    nu_ae,nu1,nu2,s,Tae,T1,m12,m21 = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate(nu_func, T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    return fs
    
    
def IM_ae_b_2M(params, ns):
    """
    nu1= pop size after pop expansion (and pop size of North Sea)
	s=proportion of the North Sea pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Baltic Sea pop
	T1= time of second epoch of migration and population growth
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
	i1= reduction in migration rate for Islands from the Baltic to North Sea 
	i2= reduction in migration rate for Islands from the North Sea to Baltic
	P= proportion of the genome made up of "islands"
    """
    nu_ae,nu1,nu2,s,Tae,T1,m12,m21,i1,i2,P = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate(nu_func, T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
# calculate teh spectrum for genomic islands	
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsi = moments.Spectrum(stsi)
    fsi.integrate([nu_ae], Tae)
    fsi = moments.Manips.split_1D_to_2D(fsi, ns[0], ns[1])
    fsi.integrate(nu_func, T1, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
	
    fs2=P*fsi+(1-P)*fs
    return fs2


######################
####   SC and SC +  ####
######################

def SC(params, ns):
    """
    nu1= pop size for North Sea
	nu2=pop size for Baltic Sea 
	T1= time of population split
	T2= time of secondary contact
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
    """
    nu1,nu2,T1,T2,m12,m21 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, m = numpy.array([[0, 0], [0, 0]]))
    fs.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    return fs


def SC_2M(params, ns):
    """
	nu1= pop size for North Sea 
	nu2= pop size for the Baltic Sea
	T1= time of population split
	T2= time of secondary contact
	i1= reduction in migration rate for Islands from the Baltic to North Sea 
	i2= reduction in migration rate for Islands from the North Sea to Baltic
	P= proportion of the genome made up of "islands"
    """
    nu1,nu2,T1,T2,m12,m21,i1,i2,P = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, m = numpy.array([[0, 0], [0, 0]]))
    fs.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
# calculate the spectrum from genomic islands
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsi = moments.Spectrum(stsi)
    fsi = moments.Manips.split_1D_to_2D(fsi, ns[0], ns[1])
    fsi.integrate([nu1, nu2], T1,  dt_fac=0.01, m = numpy.array([[0, 0], [0, 0]]))
    fsi.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fs2=P*fsi+(1-P)*fs
    return fs2


# split with growth and asymmetrical migration
def SC_ae(params, ns):
    """
    nu1= pop size following ancestral population expnasion (assumed constant in both NS and Baltic Sea)
	Tae= timing of ancestral population expansion
	T1= time of population split
	T2= time of secondary contact
	m12_0= migration rate from North Sea to Baltic
	m21_0= migration rate from Baltic Sea to North Sea
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
    """
    nu_ae,nu1,nu2,Tae,T1,T2,m12,m21 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, m = numpy.array([[0, 0], [0, 0]]))
    fs.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    return fs
    

def SC_ae_2M(params, ns):
    """
	nu1= pop size after ancestral expansion
	Tae= time of ancestral expansion
	T1= time of population split
	T2= time of secondary contact
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
	i1= reduction in migration rate for Islands from the Baltic to North Sea 
	i2= reduction in migration rate for Islands from the North Sea to Baltic
	P= proportion of the genome made up of "islands"
    """
    nu_ae,nu1,nu2,Tae,T1,T2,m12,m21,i1,i2,P = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, m = numpy.array([[0, 0], [0, 0]]))
    fs.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
# calculate the spectrum from genomic islands
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsi = moments.Spectrum(stsi)
    fsi.integrate([nu_ae], Tae)
    fsi = moments.Manips.split_1D_to_2D(fsi, ns[0], ns[1])
    fsi.integrate([nu1, nu2], T1,  dt_fac=0.01, m = numpy.array([[0, 0], [0, 0]]))
    fsi.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fs2=P*fsi+(1-P)*fs
    return fs2
    


def SC_ae_b(params, ns):
    """
    nu1= pop size after ancestral expansion (this remains constant for teh North sea population after the split)
	s=proportion of the North Sea pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Baltic Sea pop
	Tae= timing of ancestral population expansion
	T1= time of population split
	T2= time of secondary contact and start of population growth in the Baltic Sea
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
    """
    nu_ae,nu1,nu2,s,Tae,T1,T2,m12,m21 = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu1*s], T1, m = numpy.array([[0, 0], [0, 0]]))
    fs.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    return fs

def SC_ae_b_2M(params, ns):
    """
    nu1= pop size after pop expansion (and pop size of North Sea)
	s=proportion of the North Sea pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Baltic Sea pop
	T1= time of population split
	T2= time of second epoch of migration and population growth
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
	i1= reduction in migration rate for Islands from the Baltic to North Sea 
	i2= reduction in migration rate for Islands from the North Sea to Baltic
	P= proportion of the genome made up of "islands"
    """
    nu_ae,nu1,nu2,s,Tae,T1,T2,m12,m21,i1,i2,P = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu1*s], T1, m = numpy.array([[0, 0], [0, 0]]))
    fs.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
# calculate teh spectrum for genomic islands	
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsi = moments.Spectrum(stsi)
    fsi.integrate([nu_ae], Tae)
    fsi = moments.Manips.split_1D_to_2D(fsi, ns[0], ns[1])
    fsi.integrate([nu1, nu1*s], T1, m = numpy.array([[0, 0], [0, 0]]))
    fsi.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
	
    fs2=P*fsi+(1-P)*fs
    return fs2

def SC_b(params, ns):
    """
    nu1= pop size for North Sea
	s=proportion of the North Sea pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Baltic Sea pop
	T1= time of population split
	T2= time of secondary contact
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
    """
    nu1,nu2,s,T1,T2,m12,m21 = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu1*s], T1, m = numpy.array([[0, 0], [0, 0]]))
    fs.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    return fs


def SC_b_2M(params, ns):
    """
    nu1= pop size for North Sea
	s=proportion of the North Sea pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Baltic Sea pop
	T1= time of population split
	T2= time of second epoch of migration and population growth
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
	i1= reduction in migration rate for Islands from the Baltic to North Sea 
	i2= reduction in migration rate for Islands from the North Sea to Baltic
	P= proportion of the genome made up of "islands"
    """
    nu1,nu2,s,T1,T2,m12,m21,i1,i2,P = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu1*s], T1, m = numpy.array([[0, 0], [0, 0]]))
    fs.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
# calculate teh spectrum for genomic islands	
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsi = moments.Spectrum(stsi)
    fsi = moments.Manips.split_1D_to_2D(fsi, ns[0], ns[1])
    fsi.integrate([nu1, nu1*s], T1, m = numpy.array([[0, 0], [0, 0]]))
    fsi.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, m12*i2], [m21*i1, 0]]))
	
    fs2=P*fsi+(1-P)*fs
    return fs2


######################
####   AM and AM +  ####
######################

def AM(params, ns):
    """
    nu1= pop size for North Sea
	nu2=pop size for Baltic Sea 
	T1= time of population split with migration
	T2= time of speciaton
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
    """
    nu1,nu2,T1,T2,m12,m21 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, m = numpy.array([[0, m12], [0, m21]]))
    fs.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    return fs


def AM_2M(params, ns):
    """
	nu1= pop size for North Sea 
	nu2= pop size for the Baltic Sea
	T1= time of population split
	T2= time of secondary contact
	i1= reduction in migration rate for Islands from the Baltic to North Sea 
	i2= reduction in migration rate for Islands from the North Sea to Baltic
	P= proportion of the genome made up of "islands"
    """
    nu1,nu2,T1,T2,m12,m21,i1,i2,P = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, m = numpy.array([[0, m12], [m21, 0]]))
    fs.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
# calculate the spectrum from genomic islands
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsi = moments.Spectrum(stsi)
    fsi = moments.Manips.split_1D_to_2D(fsi, ns[0], ns[1])
    fsi.integrate([nu1, nu2], T1, m = numpy.array([[0, m12*i2], [m12*i1, 0]]))
    fsi.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0 ,0]]))
    fs2=P*fsi+(1-P)*fs
    return fs2


# split with growth and asymmetrical migration
def AM_ae(params, ns):
    """
    nu1= pop size following ancestral population expnasion (assumed constant in both NS and Baltic Sea)
	Tae= timing of ancestral population expansion
	T1= time of population split
	T2= time of secondary contact
	m12_0= migration rate from North Sea to Baltic
	m21_0= migration rate from Baltic Sea to North Sea
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
    """
    nu_ae,nu1,nu2,Tae,T1,T2,m12,m21 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, m = numpy.array([[0, m12], [m21, 0]]))
    fs.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    return fs
    

def AM_ae_2M(params, ns):
    """
	nu1= pop size after ancestral expansion
	Tae= time of ancestral expansion
	T1= time of population split
	T2= time of secondary contact
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
	i1= reduction in migration rate for Islands from the Baltic to North Sea 
	i2= reduction in migration rate for Islands from the North Sea to Baltic
	P= proportion of the genome made up of "islands"
    """
    nu_ae,nu1,nu2,Tae,T1,T2,m12,m21,i1,i2,P = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, m = numpy.array([[0, m12], [m21, 0]]))
    fs.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
# calculate the spectrum from genomic islands
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsi = moments.Spectrum(stsi)
    fsi.integrate([nu_ae], Tae)
    fsi = moments.Manips.split_1D_to_2D(fsi, ns[0], ns[1])
    fsi.integrate([nu1, nu2], T1, m = numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsi.integrate([nu1, nu2], T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    fs2=P*fsi+(1-P)*fs
    return fs2
    


def AM_ae_b(params, ns):
    """
    nu1= pop size after ancestral expansion (this remains constant for teh North sea population after the split)
	s=proportion of the North Sea pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Baltic Sea pop
	Tae= timing of ancestral population expansion
	T1= time of population split
	T2= time of secondary contact and start of population growth in the Baltic Sea
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
    """
    nu_ae,nu1,nu2,s,Tae,T1,T2,m12,m21 = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu1*s], T1, m = numpy.array([[0, m12], [m21, 0]]))
    fs.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
    return fs

def AM_ae_b_2M(params, ns):
    """
    nu1= pop size after pop expansion (and pop size of North Sea)
	s=proportion of the North Sea pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Baltic Sea pop
	T1= time of population split
	T2= time of second epoch of migration and population growth
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
	i1= reduction in migration rate for Islands from the Baltic to North Sea 
	i2= reduction in migration rate for Islands from the North Sea to Baltic
	P= proportion of the genome made up of "islands"
    """
    nu_ae,nu1,nu2,s,Tae,T1,T2,m12,m21,i1,i2,P = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu1*s], T1, m = numpy.array([[0, m12], [m21, 0]]))
    fs.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
# calculate teh spectrum for genomic islands	
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsi = moments.Spectrum(stsi)
    fsi.integrate([nu_ae], Tae)
    fsi = moments.Manips.split_1D_to_2D(fsi, ns[0], ns[1])
    fsi.integrate([nu1, nu1*s], T1, m = numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsi.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
	
    fs2=P*fsi+(1-P)*fs
    return fs2

def AM_b(params, ns):
    """
    nu1= pop size for North Sea
	s=proportion of the North Sea pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Baltic Sea pop
	T1= time of population split
	T2= time of secondary contact
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
    """
    nu1,nu2,s,T1,T2,m12,m21 = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu1*s], T1, m = numpy.array([[0, m12], [m21, 0]]))
    fs.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0,0], [0, 0]]))
    return fs


def AM_b_2M(params, ns):
    """
    nu1= pop size for North Sea
	s=proportion of the North Sea pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Baltic Sea pop
	T1= time of population split
	T2= time of second epoch of migration and population growth
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
	i1= reduction in migration rate for Islands from the Baltic to North Sea 
	i2= reduction in migration rate for Islands from the North Sea to Baltic
	P= proportion of the genome made up of "islands"
    """
    nu1,nu2,s,T1,T2,m12,m21,i1,i2,P = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu1*s], T1, m = numpy.array([[0, m12], [m21, 0]]))
    fs.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
# calculate teh spectrum for genomic islands	
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsi = moments.Spectrum(stsi)
    fsi = moments.Manips.split_1D_to_2D(fsi, ns[0], ns[1])
    fsi.integrate([nu1, nu1*s], T1, m = numpy.array([[0, m12*i2], [m21*i1, 0]]))
    fsi.integrate(nu_func, T2, dt_fac=0.01, m=numpy.array([[0, 0], [0, 0]]))
	
    fs2=P*fsi+(1-P)*fs
    return fs2