"""
"""
import numpy as np 

def rad2deg( angle ) : 
    # converts radians to degrees
    return angle * 180 / np.pi 

def deg2rad( angle ) :
    # converts degrees to radians
    return angle * np.pi / 180
    
class utm :
    """
    """
    def __init__(self, ref_lon, ref_height, gen_offset=True) :

        # Store internal parameters
        self.lambda0 = ref_lon
        self.h0 = ref_height
        self.a = 6378.137 * 1000 # m
        self.f = 1 / 298.257223563
        self.k0 = 0.9996
        self.E0 = 500
        self.N0 = {
            'northern':0, #m 
            'southern':10000 * 1000 # m
        }

        self.hemi = {
            'northern':1,
            'southern':-1
        }

        # Compute some preliminary values for calculations
        self.n = self.f / ( 2 - self.f )

        self.A = self.a / ( 1 + self.n ) * ( 1 + ( self.n ** 2 )/ 4 + (self.n ** 4) / 64 )

        self.alpha = [
            ( 1 / 2 ) * self.n - ( 2 / 3 ) * self.n ** 2 + ( 5/ 16 ) * self.n ** 3,
            ( 13 / 48 ) * self.n ** 2 - ( 3/5 ) * self.n ** 3,
            ( 61 / 240 ) * self.n ** 3
        ]

        self.beta = [
            (1/2) * self.n -( 2/3 ) * self.n ** 2 + ( 37 / 96 ) * self.n ** 3,
            ( 1/ 48 ) * self.n ** 2 + ( 1/15 ) * self.n ** 3,
            ( 17 / 480 ) * self.n ** 3
        ]

        self.dell = [ 
            2 * self.n - ( 2/3 ) * self.n ** 3 - 2*self.n**3,
            ( 7/3 ) * self.n ** 2 - (8/5) * self.n ** 3,
            (58/15)* self.n ** 3
        ]

        self.offset = gen_offset
        self.x0 = None
        self.y0 = None

        
        
    def llh2utm(self, llh ):
        # Using NED frame
        lat = deg2rad( llh[0] )
        lon = deg2rad( llh[1] )
        lon0 = deg2rad( self.lambda0 )
        
        if lat >= 0 :
            hemi = 'northern'
        else :
            hemi = 'southern'
        
        if len(llh) > 2 :
            h = llh[2]
            D = self.h0 - h 
        
        else :
            D = None 
        
        # Compute internmediate values
        const = (2 * np.sqrt( self.n )) / (1 + self.n )
        t = np.sinh(
            np.arctanh( np.sin(lat) ) - const * np.arctanh( const * np.sin(lat) )
        )

        dzeta = np.arctan2( t, np.cos( lon - lon0 ) )

        deta = np.arctanh( np.sin( lon - lon0 ) / ( np.sqrt( 1 + t ** 2 )) )

        sigma = 1
        tau = 0

        for j in range(3) :
            sigma += 2 * (j+1) * self.alpha[j] * np.cos( 2 * (j+1) * dzeta ) * np.cosh( 2 * (j+1) * deta )
            tau += 2 * (j+1) * self.alpha[j] * np.sin( 2 * (j+1) * dzeta) * np.sinh( 2* (j+1) * deta )
        
        # Compute the final values
        E_sum = deta
        N_sum = dzeta
        for i in range(3) :
            j = i + 1
            E_sum += self.alpha[i] * np.cos( 2*j*dzeta)*np.sinh(2*j*deta)
            N_sum += self.alpha[i] * np.sin( 2*j*dzeta)*np.cosh(2*j*deta)
        
        E = self.E0 + self.k0 * self.A * E_sum 
        N = self.N0[hemi] + self.k0 * self.A * N_sum 

        if self.offset :

            if self.x0 == None and self.y0 == None :
                self.x0 = N 
                self.y0 = E 
            
            x = N - self.x0 
            y = E - self.y0 
            
            if D != None :
                z = D 

                return (x, y, z)
            
            else:
                return (x, y)
        
        else:
            if D != None :
                return (N,E,D)
            else:
                return (N, E)
    
    def utm2llh(self, NED, hemi='northern' ) :

        # This will calculate the llh from NED given
        N = NED[0]
        E = NED[1]

        lon0 = deg2rad(self.lambda0)

        if len(NED) > 2 :
            D = NED[2]
            h = self.h0 - D 
        else :
            h = None

        # Compute intermediate values:
        zeta = (N- self.N0[hemi]) / (self.k0 * self.A)

        eta = ( E - self.E0 ) / ( self.k0 * self.A)

        dzeta = zeta 
        deta = eta
        dsigma = 1 
        dtau = 0

        for i in range(3) :
            j = i + 1

            dzeta -= self.beta[i] * np.sin( 2*j*zeta )* np.cosh( 2*j*eta)
            deta -= self.beta[i] * np.cos( 2* j * zeta) * np.sinh( 2* j * eta)
            dsigma -= 2*j*self.beta[i] * np.cos( 2*j*zeta) * np.cosh( 2*j*eta )
            dtau += 2*j*self.beta[i] * np.sin( 2* j* zeta) * np.sinh(2 * j * eta)
        
        chi = np.arcsin( np.sin(dzeta) / np.cosh( deta ) )

        # Compute the final values
        lat = chi 

        for i in range(3) :
            j = i + 1
            lat+= self.dell[i] * np.sin( 2 * j* chi)
        
        lon = lon0 + np.arctan( np.sinh(deta) / np.cos( dzeta ) )

        # Reconvert to degrees
        lat = rad2deg( lat )
        lon = rad2deg( lon )

        if h != None :
            return ( lat, lon, h )
        else:
            return ( lat, lon )

    
