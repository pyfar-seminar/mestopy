import warnings
import numpy as np
from pyfar import Signal


# function to perform a deconvolution of two signals
def deconv(x_meas, x_ref, measurement_chain=None):
    # Check if both inputs are type Signal
    if not isinstance(x_meas, Signal):
        raise TypeError('Input data has to be of type: Signal.')
    if not isinstance(x_ref, Signal):
        raise TypeError('Input data has to be of type: Signal.')

    # Transform both signals to time domain
    x_meas.domain = 'time'
    x_ref.domain = 'time'

    # Check if both signals have the same sampling rate
    if x_meas.sampling_rate != x_ref.sampling_rate:
        raise ValueError("The two signals have different sampling rates!")
    # Check if both signals have the same fft norm, if not: warn
    if x_meas.fft_norm != x_ref.fft_norm:
        warnings.warn("The two signals have different fft_norms.")

    # Check if both signals have the same length,
    # if not: bring them to the same length
    if x_meas.n_samples > x_ref.n_samples:
        # Add Zeros to x_ref
        x_ref = Signal(np.concatenate((x_ref.time,
                                       np.zeros(x_meas.n_samples -
                                                x_ref.n_samples))),
                       x_ref.sampling_rate, fft_norm=x_ref.fft_norm,
                       dtype=x_ref.dtype, comment=x_ref.comment)

    if x_meas.n_samples < x_ref.n_samples:
        # Add Zeros to x_meas
        x_meas = Signal(np.concatenate((x_meas.time,
                                        np.zeros(x_ref.n_samples -
                                                 x_meas.n_samples))),
                        x_meas.sampling_rate,
                        fft_norm=x_meas.fft_norm,
                        dtype=x_meas.dtype,
                        comment=x_meas.comment)

    # Transform both signals to frequency domain and
    # divide to get the transfer funktion
    x_meas.domain = 'freq'
    x_ref.domain = 'freq'
    H = x_meas/x_ref

    # Check if the signals have any comments,
    # if yes: concatenate the comments for the result
    if not (x_meas.comment is None and x_ref.comment is None):
        H.comment = ('IR calculated with deconvolution: [1]' +
                     x_meas.comment + '; [2]' + x_ref.comment)
    else:
        H.comment = "IR calculated with deconvolution"

    # Transform back to time domain and return the impulse resonse
    # FFT-Norm auf none setzen
    H.domain = 'time'
    return H


# function to perform complete measurement and save as pyfar.Signal
def take(x_excitation, measurement_chain=None):
    # play signal
    # record signal
    # deconv
    # + cal resp
    """TODO: to be defined. """
    return 0  # pyfar-Signal of Impulse Response


# Class to generate ref-Objects, that can bei part of the MeasurementChain
class DeviceObj(object):

    def __init__(self, data, sens=1, name=''):
        if not isinstance(data, Signal):
            raise AttributeError("Useless DeviceObj,"
                                 " data or sens must be defined!")
        self.sens = sens
        self.name = name
        self.data = data

    @property
    def freq(self):
        self.data.domain = 'freq'
        return self.data * self.sens

    @property
    def device_name(self):
        return self.name

    @device_name.setter
    def device_name(self, new_name):
        self.name = new_name

    def __repr__(self):
        """String representation of DeviceObj class.
        """
        repr_string = (
            f"{self.name} defined by {self.data.n_bins} freq-bins, "
            f"sensitivity={self.sens}\n")
        return repr_string


# Class for MeasurementChain as frame for RefObjs and calibration
class MeasurementChain(object):

    def __init__(self,
                 sampling_rate,
                 sound_device=None,
                 devices=[],
                 comment=None):
        self.sampling_rate = sampling_rate
        self.sound_device = sound_device
        self.devices = devices
        self.comment = comment

    def add_device(self,
                   device_data=None,
                   sens=1,
                   device_name=''):
        if device_data is None:
            device_data = Signal(np.full(int(self.sampling_rate/2), 1.0),
                                 self.sampling_rate,
                                 domain='freq')
        # check if ref_signal is a pyfar.Signal, if not raise Error
        if not isinstance(device_data, Signal):
            raise TypeError('Input data must be of type: Signal.')
        # check if there are no devices in measurement chain
        if self.devices == []:
            # add ref-measurement to chain
            device_data.domain = 'freq'
            new_device = DeviceObj(device_data,
                                   sens,
                                   device_name)
            self.devices.append(new_device)
        else:
            # check if n_bins of all devices is the same
            if not self.devices[0].data.n_bins == device_data.n_bins:
                raise ValueError("ref_signal has wrong n_bins")
            # check if sampling_rate of new device and MeasurementChain
            # is the same
            if not self.sampling_rate == device_data.sampling_rate:
                raise ValueError("ref_signal has wrong samping_rate")
            # add device to chain
            new_device = DeviceObj(device_data,
                                   sens,
                                   device_name)
            self.devices.append(new_device)

    def list_devices(self):
        # list all ref-objects in chain
        device_names = []
        for dev in self.devices:
            name = dev.device_name
            device_names.append(name)
        return device_names

    def remove_device(self, num):
        # remove ref-object in chain position num
        if isinstance(num, int):
            front = self.devices[:num]
            back = self.devices[num+1:]
            for i in back:
                front.append(i)
            self.devices = front
        # remove ref-object in chain by name
        elif isinstance(num, str):
            i = 0
            for dev in self.devices:
                if dev.name == num:
                    self.remove_device(i)
                    return
                i = i + 1
            raise ValueError(f"device {num} not found")
        else:
            raise TypeError("device to remove must be int or str")

    # reset complete ref-object-list
    def reset_devices(self):
        self.devices = []

    # get the freq-response of specific device in measurement chain
    def device_freq(self, num):
        if isinstance(num, int):
            dev = self.devices[num].data * self.devices[num].sens
            dev.domain = 'freq'
            return dev
        elif isinstance(num, str):
            i = 0
            for dev in self.devices:
                if dev.name == num:
                    return self.device_freq(i)
                i = i + 1
            raise ValueError(f"device {num} not found")
        else:
            raise TypeError("device to remove must be int or str")

    # get the freq-response of whole measurement chain as pyfar.Signal
    def freq(self):
        if self.devices != []:
            resp = Signal(np.ones(int(self.devices[0].data.n_bins)),
                          self.sampling_rate,
                          domain='freq',
                          fft_norm=self.devices[0].data.fft_norm,
                          dtype=self.devices[0].data.dtype)
            for dev in self.devices:
                resp = resp * dev.data * dev.sens
        else:
            resp = Signal(np.ones(self.sampling_rate), self.sampling_rate)
        return resp

    def __repr__(self):
        """String representation of MeasurementChain class.
        """
        repr_string = (
            f"measurement chain with {len(self.devices)} devices "
            f"@ {self.sampling_rate} Hz sampling rate.\n")
        i = 1
        for dev in self.devices:
            repr_string = f"{repr_string}# {i}: {dev}"
            i += 1
        return repr_string
