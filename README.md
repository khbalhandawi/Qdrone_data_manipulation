# Qdrone_data_manipulation

Python file is used to import raw QDrone data "optitrack_Data.m" and perform the following operations:

- Load data from .mat file
- Parse data for differentiation operation
- Differentiate parsed data
- Unparse the differentiated position data to obtain velocity signal in sync with IMU measurements
- Inroduce an artificial delay in the position signal given by the following options (line 407):
	- delay = 20; #samples
- Downsample raw position input using the follow options (line 414):
	- downsample_frequency_units = 40   -----> downsample 100 Hz raw signal to 25 Hz
	- downsample_frequency_units = 200  -----> downsample 100 Hz raw signal to 5 Hz
	- downsample_frequency_units = 1000 -----> downsample 100 Hz raw signal to 1 Hz