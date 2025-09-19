import atlite
import os
import yaml
import argparse


### --------- Input data and configuration ------------- ###

dirname = os.getcwd() 
with open(os.path.join("configs/config.yaml"), "r", encoding="utf-8") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

weather_year = config["weather_year"]
weather_data_path = os.path.abspath(config["weather_data_path"])
bounds = config["weather_data_area_bounds"]
country_code = config["country_code"]

# Initialize parser for command line arguments and define arguments
parser = argparse.ArgumentParser()
parser.add_argument("--method",default="manual", help="method to run the script, e.g., snakemake or manual")
parser.add_argument("--weather_year", default=weather_year, help="weather year for the energy profiles")
args = parser.parse_args()

# If running via Snakemake, use the region name and folder name from command line arguments
if args.method == "snakemake":
    weather_year = args.weather_year
    print(f'\nWeather data preparation')
    print(f"Running via snakemake - measures: weather_year={weather_year}")
else:
    print('\nWeather data preparation')
    print(f"Running manually - measures: weather_year={weather_year}")

#----------- Prepare weather data cutouts from ERA5 via Copernicus API ------------- ###
print("Processing weather year: ", weather_year)

t_start_1 = f"{weather_year}-01-01"
t_end_1 = f"{weather_year}-12-31"

# Define cutouts
cutout_file_path = os.path.join(weather_data_path, f"{country_code}_{weather_year}")
cutout = atlite.Cutout(cutout_file_path, module="era5", x=slice(bounds['west'], bounds['east']), y=slice(bounds['south'], bounds['north']), time=slice(t_start_1, t_end_1))

# Connect to API and preprare/download cutouts
print("Preparing cutout from API and downloading data ...")
cutout.prepare()