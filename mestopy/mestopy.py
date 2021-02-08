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
def take(x_excitation, meas_chain):
    # play signal
    # record signal
    # deconv
    # + cal resp
    """TODO: to be defined. """
    return 0  # pyfar-Signal of Impulse Response


# Class to generate ref-Objects, that can bei part of the MeasurementChain
class RefObj(object):

    def __init__(self, data, calibration=1, name='', inverted=False):
        self.data = data
        self.calibration = calibration
        self.name = name
        self.inverted = inverted
        """TODO: to be defined. """

    @property
    def freq(self):
        self.data.domain = 'freq'
        return self.data

    @freq.setter
    def freq(self, ref_signal, calibration=1, device_name='', inverted=False):
        ref_signal.domain = "freq"
        self.data = ref_signal
        self.name = device_name
        self.inverted = inverted

    @property
    def device_name(self):
        return self.name

    @device_name.setter
    def device_name(self, new_name):
        self.name = new_name

    @property
    def is_inverted(self):
        return self.inverted

    @is_inverted.setter
    def invert(self):
        new_data = self.ref_data.freq() * (-1)
        self.data = Signal(new_data,
                           self.data.sampling_rate(),
                           domain='freq',
                           fft_norm=self.data.fft_norm(),
                           dtype=self.data.dtype())


# Class for MeasurementChain as frame for RefObjs and calibration
class MeasurementChain(object):

    def __init__(self, sampling_rate, sound_device, refs=[], comment=None):
        self.sampling_rate = sampling_rate
        self.sound_device = sound_device
        self.refs = refs
        self.comment = comment
        """TODO: to be defined. """

    def add_ref(self,
                ref_signal,
                calibration=1,
                device_name='',
                inverted=False):
        # check if ref_signal is a pyfar.Signal, if not raise Error
        if not isinstance(ref_signal, Signal):
            raise TypeError('Input data has to be of type: Signal.')
        if self.refs == []:
            # add ref-measurement to chain
            ref_signal.domain = 'freq'
            new_ref = RefObj(ref_signal,
                             calibration,
                             device_name,
                             inverted=inverted)
            self.refs.append(new_ref)
        else:
            # check if n_bins of all ref_signals is the same
            if (self.refs[0].data.n_bins == ref_signal.n_bins and
                    (self.refs[0].data.sampling_rate ==
                     ref_signal.sampling_rate)):
                # add ref-measurement to chain
                new_ref = RefObj(ref_signal,
                                 calibration,
                                 device_name,
                                 inverted=inverted)
                self.refs.append(new_ref)
            else:
                raise ValueError("ref_signal has wrong samping_rate or n_bins")

    def ls_ref(self):
        # list all ref-objects in chain
        ref_names = []
        for ob in self.refs:
            name = ob.device_name
            ref_names.append(name)
        return ref_names

    def rm_ref(self, num):
        # remove ref-object in chain position num
        if isinstance(num, int):
            front = self.refs[:num]
            back = self.refs[num+1:]
            for i in back:
                front.append(i)
            self.refs = front
        else:
            raise TypeError("ref-object to remove must be int")

    # reset complete ref-object-list
    def reset_ref(self):
        self.refs = []

    # get the freq-response of specific device in measurement chain
    def get_ref(self, num):
        if isinstance(num, int):
            resp = self.refs[num].data
            resp.domain = 'freq'
            return resp
        else:
            raise TypeError("ref-object to remove must be int")

    # get the freq-response of whole measurement chain as pyfar.Signal
    def get_refs(self):
        if self.refs != []:
            resp = Signal(np.ones(int(self.refs[0].data.n_bins)),
                          self.sampling_rate,
                          domain='freq',
                          fft_norm=self.refs[0].data.fft_norm,
                          dtype=self.refs[0].data.dtype)
            for ref in self.refs:
                resp = resp * ref.data
        else:
            resp = Signal(np.ones(self.sampling_rate),
                          self.sampling_rate)
        return resp
