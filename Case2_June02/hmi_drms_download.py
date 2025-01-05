import os
import drms


#email_address = os.environ["adithya1@atulbhats.com"]  
client = drms.Client(email='adithya1@atulbhats.com') 
out_dir='/Analysis/Projects_Data/Flare_Data/June02_Flare_Data1/hmi2'
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

export_request = client.export('hmi.M_45s[2024.06.02_TAI/1d@45s]{Magnetogram}', method='url', protocol='fits')
export_request.wait()
export_request.download(out_dir, index=0) 