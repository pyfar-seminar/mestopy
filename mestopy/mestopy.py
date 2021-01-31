import warnings
import numpy as np
from pyfar import Signal


# function to perform a deconvolution of two signals
def deconv(x_meas, x_ref):
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
def take(x_excitation, meas_chain):
    # play signal
    # record signal
    # deconv
    # + cal resp

    return 0  # pyfar-Signal of Impulse Response


# Class to generate ref-Objects, that can bei part of the MeasurementChain
class RefObj(object, Signal):

    def __init__(
            self,
            ref_data,
            calibration=1,
            device='',
            inverted=False):
        """TODO: to be defined. """

    @property
    def freq(self):
        self.ref_data.domain = 'freq'
        return self.ref_data._data

    @freq.setter
    def freq(self, ref_signal, calibration=1, device_name='', inverted=False):
        ref_signal.domain = "freq"
        self._ref_data = ref_signal
        self._device = device_name
        self._inverted = inverted

    @property
    def device_name(self):
        return self._device

    @device_name.setter
    def device_name(self, new_name):
        self._device = new_name

    @property
    def inverted(self):
        return self._inverted

    @inverted.setter
    def invert(self):
        new_data = self.ref_data.freq() * (-1)
        self._ref_data = Signal(new_data,
                                self.ref_data.sampling_rate(),
                                self.ref_data.sampling_rate())


# Class for MeasurementChain as frame for RefObjs and calibration
class MeasurementChain(object, Signal, RefObj):

    def __init__(self,
                 sampling_rate,
                 sound_device,
                 refs=[],
                 comment=None):
        """TODO: to be defined. """

    def add_ref(self, ref_signal, calibration, device_name, inverted):
        # add ref-measurement to chain
        new_ref = RefObj()
        new_ref.freq(ref_signal, calibration, device_name, inverted)
        self.refs.append(new_ref)

    def ls_ref(self):
        # list all ref-objects in chain
        ref_names = []
        for ob in self.calibration_responses:
            name = ob.device_name()
            ref_names.append(name)
        return ref_names

    def rm_ref(self, num):
        # remove ref-object in chain position num
        if isinstance(num, int):
            pass
        # remove ref-object in chain by name
        elif isinstance(num, str):
            pass
        else:
            raise TypeError("ref-object to remove must be int or str!")
