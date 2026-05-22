import numpy as np
from astropy.time import Time
import astropy.units as u

def find_cotemporal_peaks(suit_times, helios_times, tolerance_sec=90):
    """
    suit_times : astropy Time array of SUIT peak times
    helios_times : astropy Time array of Helios peak times
    tolerance_sec : time window for matching
    """
    
    matches = []
    tolerance = tolerance_sec * u.second
    
    for st in suit_times:
        dt = np.abs(helios_times - st)
        min_idx = np.argmin(dt)
        
        if dt[min_idx] < tolerance:
            matches.append((st, helios_times[min_idx]))
    
    return matches

# Suppose you already extracted peak indices

helios=np.loadtxt("helios_peaks.csv", delimiter=',', skiprows=1, dtype='str')
suit=np.loadtxt("suit_diff_peaks.csv", delimiter=',', skiprows=1, dtype='str')
suit_time_array = np.array(suit[:,0], dtype='datetime64[s]')
helios_time_array = np.array(helios[:,0], dtype='datetime64[s]')
suit_counts=np.array(suit[:,1],dtype=float)
helios_counts=np.array(helios[:,1],dtype=float)

suit_peak_times = Time(suit_time_array)
helios_peak_times = Time(helios_time_array)

matches = find_cotemporal_peaks(suit_peak_times,
                                helios_peak_times,
                                tolerance_sec=97)

print("Number of co-temporal peaks:", len(matches))
# print(matches)
def get_corresponding_counts(match_times, helios_times,suit_counts,helios_counts):
    hcounts = []
    scounts = []
    #print(match_times[0],helios_times[0],suit_peak_times[0])
    for t in match_times:
        ht = np.abs(helios_times - t[0])
        hidx = np.argmin(ht)
        hcounts.append(helios_counts[hidx])
        st = suit_peak_times[np.argmin(np.abs(suit_peak_times - t[1]))]
        sidx = np.argmin(np.abs(suit_peak_times - st))
        scounts.append(suit_counts[sidx])  # assuming second column is counts
    
    

    return np.array(hcounts), np.array(scounts)

hcounts, scounts = get_corresponding_counts(matches, helios_peak_times,suit_counts,helios_counts)
print(hcounts,scounts)
# print(type(matches))

suit_dt_arr = np.array([(t1.iso) for t1, t2 in matches])
np.savetxt('suit_helios_time_count_c07.csv',np.c_[suit_dt_arr,scounts,hcounts],header='date_time,suit_counts,helios_counts',comments='',delimiter=',',fmt='%s')