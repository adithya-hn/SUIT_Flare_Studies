
from sunpy.net import Fido
from sunpy.net import attrs as a

event_type = "FL"
tstart = "2024-11-13 00:10:00"
tend = "2024-11-13 20:00:00"
result = Fido.search(a.Time(tstart, tend),
                     a.hek.EventType(event_type),
                     a.hek.FL.GOESCls > "M0.0",
                     a.hek.OBS.Observatory == "GOES")
#print(result.show("hpc_bbox", "refs"))
print(result.show("event_starttime", "event_peaktime", "event_endtime", "hpc_bbox", "fl_goescls", "obs_observatory"))

'''from sunpy.net import hek
client = hek.HEKClient()
results = client.query(hek.attrs.Time("2025-09-01", "2025-09-24"),
                       hek.attrs.EventType("FL"))'''