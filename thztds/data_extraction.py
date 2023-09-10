"""
# data_extraction
This is used to extract different data values from THz-TDS (Terahertz time-domain spectrocopy)

The default units are:
* time [pico-seconds] 
* frequency [Tera-hertz]
* distane [micro-meters]
"""
import pandas as pd
import numpy as np
import json
import thztds.spectrum as spectrum
import scipy.constants as constant

light_speed = constant.speed_of_light * 1e-6  # µm/ps


def import_data(file) -> tuple:
    """
    Import the data from the file giving the time and amplitude as function of time
    ### Input
    * `file` = file path
    ### Output
    * `(time, amplitude)` = the time and amplitude as lists of data
    """
    # Extract data from file
    data = pd.read_csv(file, header=None, sep="\t", decimal=",")

    # gives time and amplitude
    time, amplitude = data.iloc[:, 0], data.iloc[:, 1]

    # Transform into lists
    time = time.to_list()
    amplitude = amplitude.to_list()

    return time, amplitude


def spectral_decomposition(amplitude, time_step) -> tuple:
    """
    Give the frequency and the complex amplitude as function of frequency
    ### Input
    `time_step` = The time difference between samples
    `amplitude` = The amplitude data as function of time

    ### Output
    `(frequency,amplitude)`
    """
    amplitude = np.fft.fft(amplitude)
    n = amplitude.size  # The number of elemnets in the data of amplitude
    frequency = np.fft.fftfreq(n, d=time_step)

    # Remove the negative frequencies
    pos = len(amplitude) // 2
    frequency = frequency[:pos]
    amplitude = amplitude[:pos]

    return frequency, amplitude


def calculate_transmission_spectrum(
    file_sample, file_reference, time_step=None
) -> dict:
    """
    ### Input
    `file_sample` = path to the file with the sample data
    `file_reference` = path to the file with reference data
    `time_step` = time step value
    """
    # Import the time and amplitude data
    time_sample, amplitude_sample_time = import_data(file_sample)
    time_reference, amplitude_reference_time = import_data(file_reference)

    # Convert the infromation to frequency and amplitude
    frequency_sample, amplitude_sample_frequency = spectral_decomposition(
        amplitude_sample_time, time_step
    )
    frequency_reference, amplitude_reference_frequency = spectral_decomposition(
        amplitude_reference_time, time_step
    )

    # Calculate the transmission from the two data
    transmission = amplitude_sample_frequency / amplitude_reference_frequency
    phase_transmission = np.unwrap(np.angle(transmission))

    return {
        "time_sample": time_sample,
        "amplitude_sample_time": amplitude_sample_time,
        "time_reference": time_reference,
        "amplitude_reference_time": amplitude_reference_time,
        "frequency_sample": frequency_sample,
        "amplitude_sample_frequency": amplitude_sample_frequency,
        "frequency_reference": frequency_reference,
        "amplitude_reference_frequency": amplitude_reference_frequency,
        "transmission": transmission,
        "phase_transmission": -(phase_transmission - phase_transmission[0]),
    }


def create_data_file(file, *time_step):
    """
    create a data file
    ### Input
    * `file` = file path with extension
    * `*time_step` = the time step between consecutive values
    """
    time, amplitude = import_data(file)
    if not (time_step):
        time_step = time[1] - time[0]
    frequency, amplitude_frequency = spectral_decomposition(amplitude, time_step)

    time, amplitude = pd.DataFrame(time), pd.DataFrame(amplitude)
    data = pd.concat([time, amplitude, frequency, amplitude_frequency])

    data.to_csv(file, sep="\t", decimal=".", header=True)


def extract(json_file: str, output_file: str):
    file = open(json_file)
    data = json.load(file)

    data_info = data["material"]

    name = data_info["name"]
    thickness = data_info["thickness"]
    sample = data_info["sample"]
    reference = data_info["reference"]
    time_step = data_info["time_step"]
    
    initial_refractive_index = (
        data_info["initial_refractive_index"][0]
        + 1j * data_info["initial_refractive_index"][1]
    )

    file.close()

    # Import the time and amplitude data
    time_sample, amplitude_time_sample = import_data(sample)
    time_reference, amplitude_time_reference = import_data(reference)


    csv_file_data_time = pd.DataFrame(
        {
            "time_sample": time_sample,
            "amplitude_time_sample": amplitude_time_sample,
            "time_reference": time_reference,
            "amplitude_time_reference": amplitude_time_reference,
        }
    )

    # Calculate the frequency and amplitude data
    frequency_sample, amplitude_frequency_sample = spectral_decomposition(
        amplitude_time_sample, time_step
    )
    frequency_reference, amplitude_frequency_reference = spectral_decomposition(
        amplitude_time_reference, time_step
    )

    # Calculate the transmission from the sample and reference data
    transmission = amplitude_frequency_sample / amplitude_frequency_reference
    phase_transmission = np.unwrap(np.angle(transmission))

    csv_file_data_frequency = pd.DataFrame(
        {
            "frequency_sample": frequency_sample,
            "amplitude_frequency_sample_re": np.real(amplitude_frequency_sample),
            "amplitude_frequency_sample_im": np.imag(amplitude_frequency_sample),
            "frequency_reference": frequency_reference,
            "amplitude_frequency_reference_re": np.real(amplitude_frequency_reference),
            "amplitude_frequency_reference_im": np.imag(amplitude_frequency_reference),
            "transmission_re": np.real(transmission),
            "transmission_im": np.imag(transmission),
            "phase_transmission": -(phase_transmission),
        }
    )

    refractive_index = spectrum.Medium(name, thickness).refractive_index(
        transmission,
        list(csv_file_data_frequency.loc[:, "phase_transmission"]),
        2 * np.pi * frequency_sample,
        initial_refractive_index,
        5,
    )

    absorption_coefficient = []

    for i, j in zip(refractive_index, frequency_sample):
        absorption_coefficient.append(4 * np.pi * np.imag(i) * j * 1e4 / light_speed)

    csv_file_data_refractive_index = pd.DataFrame(
        {
            "refractive_index_re": np.real(refractive_index),
            "refractive_index_im": np.imag(refractive_index),
            "absorption_coefficient": absorption_coefficient,
        }
    )

    n_c = np.array(refractive_index, dtype=complex)
    f = np.array(frequency_sample, dtype=float)
    n = np.real(n_c)
    k = np.imag(n_c)
    dielectric_constant = (n**2 - k**2) + 1j * 2 * n * k
    conductivity = dielectric_constant*f/(2j)
    csv_file_data_conductivity = pd.DataFrame(
        {
            "dielectric_constant_re": np.real(dielectric_constant),
            "dielectric_constant_im": np.imag(dielectric_constant),
            "conductivity_re": np.real(conductivity),
            "conductivity_im": np.imag(conductivity),
        }
    )

    csv_file_data_time = pd.DataFrame(csv_file_data_time)
    csv_file_data_frequency = pd.DataFrame(csv_file_data_frequency)

    # Concatenate all the data information
    csv_file_data = pd.concat(
        [
            csv_file_data_time,
            csv_file_data_frequency,
            csv_file_data_refractive_index,
            csv_file_data_conductivity,
        ],
        axis=1,
    )

    # Generate de csv file
    csv_file_data.to_csv(output_file)
