import numpy as np
import warnings
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from functools import partial

def lacey_mixing_curve_extended(time, k, tau, m0, m_plateau=1.0):
    """
    ## Original DEMToolbox function description from Jack Grogan

    ### Curve for the Lacey mixing model.

    The Lacey mixing model is a exponential model that describes the 
    mixing of a binary system. The model is classically defined as [1] 
    [2]:

    .. math:: 
    
        M(t) = 1 - (1 - M_0) e^{-kt}

    where :math:`M(t)` is the Lacey mixing index at time :math:`t`, 
    :math:`M_0` is the initial Lacey Mixing index, :math:`k` is the 
    mixing rate constant, and :math:`t` is the time.

    An extension of the Lacey mixing model to include a time constant
    :math:`\\tau` to allow for delayed onset of exponential mixing 
    mixing behavior was proposed by Ratnayake et al. [2]:

    .. math::

        M(t) = max(1 - (1 - M_0) e^{-k(t - \\tau)})

    This mixing model is used to fit the mixing data as with no delayed
    onset of exponential mixing Ratnayake et al.'s model reduces to the
    classical Lacey mixing model.

    ### References

    [1] Lacey PM. Developments in the theory of particle mixing. 
        Journal of applied chemistry. 1954 May;4(5):257-68.

    [2] Ratnayake P, Chandratilleke R, Bao J, Shen Y. A soft-sensor 
        approach to mixing rate determination in powder mixers. Powder 
        Technology. 2018 Aug 1;336:493-505.

    ## MEng-RP additions

    In my MEng-RP a singular value of resolution was used to sample multiple 
    fill heights, as such  the bin volumes vary across runs. This may mean that 
    :math:`M_0 != 0` and varies across runs, and also may mean that 
    :math:`M_{final} != 1`. The above proposed model already uses :math:`M_0` 
    as a parameter, and as such no changes are made in that respect. However, the 
    final value of the Lacey mixing index is taken to be 1 in these models. For 
    the purpose of time-series analysis in my MEng-RP, the term :math:`M_{plateau}` 
    is used to represent the value of Lacey mixing index when the system cannot 
    become more mixed, where :math:`0<M_{plateau}<=1.` As such the Lacey mixing 
    model and extended model are expressed as:

    .. math::
    
        M(t) = M_{plateau} - (M_{plateau} - M_0) e^{-kt}.

    and,

    .. math::
        M(t) = max(M_{plateau} - (M_{plateau} - M_0) e^{-k(t - \\tau)})
    
    respectively.
    
    Parameters
    ----------
    time : array-like
        The time data for the lacey mixing index.
    k : float
        The rate constant with units of inverse time.
    tau : float
        The time constant. Allows for delayed onset of exponential 
        mixing.
    m0 : float
        The initial value of the Lacey mixing index.
    m_plateau : float, optional
        The plateau value of the Lacey mixing index, by default 1.0.

    Returns
    -------
    array-like
        The predicted values for the Lacey mixing curve.

    Raises
    ------
    ValueError
        If time is not an array-like object.
    ValueError
        If k is not an integer or float.
    ValueError
        If tau is not an integer or float.
    ValueError
        If m0 is not an integer or float.
    """
    if not isinstance(time, np.ndarray):
        time = np.array(time)
        raise ValueError("time must be an array-like object")
    
    if not isinstance(k, (int, float)):
        raise ValueError("k must be an integer or float")
    
    if not isinstance(tau, (int, float)):
        raise ValueError("tau must be an integer or float")
    
    if not isinstance(m0, (int, float)):
        raise ValueError("m0 must be an integer or float")
    
    if not isinstance(m_plateau, (int, float)) or not (0 < m_plateau <= 1):
        raise ValueError("m_plateau must be an integer or float between 0 (exclusive) "
        "and 1 (inclusive)")
    
    return [max((m_plateau - (m_plateau - m0) * np.exp(-k*(t - tau))), m0) for t in time]

def lacey_mixing_curve_fit_extended(time, m, window_size=20, t0=0, tend=None, k0=0.1, free_plateau=False):
    """
    ## Original DEMToolbox function description from Jack Grogan

    ### Fit the Lacey mixing curve to the data.

    The Lacey mixing model is a exponential model that describes the 
    mixing of a binary system. The model is classically defined as [1] 
    [2]:

    .. math:: 
    
        M(t) = 1 - (1 - M_0) e^{-kt}

    where :math:`M(t)` is the Lacey mixing index at time :math:`t`, 
    :math:`M_0` is the initial Lacey Mixing index, :math:`k` is the 
    mixing rate constant, and :math:`t` is the time.

    An extension of the Lacey mixing model to include a time constant
    :math:`\\tau` to allow for delayed onset of exponential mixing 
    mixing behavior was proposed by Ratnayake et al. [2]:

    .. math::

        M(t) = max(1 - (1 - M_0) e^{-k(t - \\tau)})

    This mixing model is used to fit the mixing data as with no delayed
    onset of exponential mixing Ratnayake et al.'s model reduces to the
    classical Lacey mixing model.

    ### References

    [1] Lacey PM. Developments in the theory of particle mixing. 
        Journal of applied chemistry. 1954 May;4(5):257-68.

    [2] Ratnayake P, Chandratilleke R, Bao J, Shen Y. A soft-sensor 
        approach to mixing rate determination in powder mixers. Powder 
        Technology. 2018 Aug 1;336:493-505.

    ## MEng-RP additions

    In my MEng-RP a singular value of resolution was used to sample 
    multiple fill heights, as such  the bin volumes vary across runs. 
    This may mean that :math:`M_0 != 0` and varies across runs, and also 
    may mean that :math:`M_{final} != 1`. The above proposed models already 
    use :math:`M_0` as a parameter, and as such no changes are made in that 
    respect. However, the final value of the Lacey mixing index is taken to 
    be 1 in these models. For the purpose of time-series analysis in my MEng-RP, 
    the term :math:`M_{plateau}` is used to represent the value of Lacey 
    mixing index when the system cannot become more mixed, where 
    :math:`0<M_{plateau}<=1.` As such the Lacey mixing model and the extended 
    model are expressed as:

    .. math::
    
        M(t) = M_{plateau} - (M_{plateau} - M_0) e^{-kt},

    and
    .. math::
        M(t) = max(M_{plateau} - (M_{plateau} - M_0) e^{-k(t - \\tau)})

    respectively.

    Due to the nature of sampling, :math:`M` continues to vary slightly, even after 
    a "stable" degree of mixing has been reached. :math:`M_{plateau}` was calculated 
    as the mean of the final `window_size` values of LMI in the timeseries.
    In my MEng-RP simulations, this was 20, as this was the point when vibration
    stopped.

    ### Mixing Time, :math:`t_{95}`

    The mixing time of a system can be characterised using the time taken to reach 95%
    of the maximum degree of mixing, :math:`M_{95}`, and is denoted as :math:`t_{95}`. 
    :math:`M_{95}` can be calculated using the parameters from the fit as:

    .. math::
        M_{95} = 0.95 * M_{plateau}

    The :math:`t_{95}` of the input data is the time at which the first value of the 
    data set reaches or exceeds :math:`M_{95}`, if this is not reached then None is returned. 
    The :math:`t_{95}` of the fitted curve is calculated by rearranging the model of t, for a
    value of M_{95}:

    .. math::
        t_{95} = \\tau - \\frac{1}{k} \\ln \\left( \\frac{0.05 * M_{plateau}}{M_{plateau} - M_0} \\right)


    Parameters
    ----------
    time : array-like
        The time data for the lacey mixing index.
    m : array-like
        The lacey mixing index data.
    window_size : int, optional
        The number of consecutive values to consider when determining if the 
        system has reached the plateau value of the Lacey mixing index, by default 20.
    t0 : int or float
        The time at which mixing begins by default 0.
    tend : int or float, optional
        The time at which mixing ends, by default None. If None, then
        all the time data from the start time will be used.
    k0 : float, optional
        The initial guess for the rate constant k, by default 0.1.
    free_plateau : bool, optional
        If True, m_plateau is treated as a free parameter during fitting, using
        the empirical value as an initial guess. Defaults to False.

    Returns
    -------
    popt : array-like
        The optimal values for the parameters k and tau (and m_plateau if
        free_plateau=True).
    pcov : 2D array
        The estimated covariance of popt as returned by curve_fit.
    time_mixing : array-like
        The time data for the mixing period used in the fit.
    m_mixing : array-like
        The lacey mixing index data for the mixing period used in the fit.
    m_fit : array-like
        The predicted lacey mixing index values for the mixing period calculated 
        using the optimal parameters.
    m_plateau : float
        The calculated plateau value of the Lacey mixing index (empirical
        average if free_plateau=False, optimised value if free_plateau=True).
    r2 : float
        The coefficient of determination of the fit.
    t95_data : float or None
        The first time in the actual data where the Lacey mixing index reaches 
        or exceeds 95% of m_plateau. None if not reached.
    t95_fit : float or None
        The interpolated time in the fitted curve where the Lacey mixing index reaches 
        95% of m_plateau.

    Raises
    ------
    ValueError
        If time is not an array-like object.
    ValueError
        If m is not an array-like object.
    ValueError
        If window_size is not an integer > 0.
    ValueError
        If std_threshold is not a float or integer > 0.
    ValueError
        If t0 is not an integer or float.
    ValueError
        If tend is not an integer or float.
    ValueError
        If time and m are not the same length.
    ValueError
        If k0 is not an integer or float.
    ValueError
        If fewer than 3 valid data points remain after preprocessing.

    Warnings
    --------
    UserWarning
        If fewer data points than `window_size`:
        Mean of all available LMI values used for m_plateau.
    UserWarning
        If LMI value in final `window_size` sample falls outside ± 1 standard deviation of mean:
        User message raised.
    """
    if not isinstance(time, np.ndarray):
        time = np.array(time)
        raise ValueError("time must be an array-like object")
    
    if not isinstance(m, np.ndarray):
        m = np.array(m)
        raise ValueError("m must be an array-like object")
    
    if not isinstance(window_size, int) or window_size <= 0:
        raise ValueError("window_size must be an integer > 0")

    if not isinstance(t0, (int, float)):
        raise ValueError("t0 must be an integer or float")
    
    if tend is None:
        tend = time[-1]
    elif not isinstance(tend, (int, float)):
        raise ValueError("tend must be an integer or float")
    
    if len(time) != len(m):
        raise ValueError("time and m must be the same length")
    
    if not isinstance(k0, (int, float)):
        raise ValueError("k0 must be an integer or float")
    
    mixing_indices = (time >= t0) & (time <= tend)
    time_mixing = time[mixing_indices]
    m_mixing = m[mixing_indices]

    # Preprocessing: drop NaNs
    valid_mask = ~np.isnan(m_mixing) & ~np.isnan(time_mixing)
    time_mixing = time_mixing[valid_mask]
    m_mixing = m_mixing[valid_mask]

    # Sort by time
    sort_idx = np.argsort(time_mixing)
    time_mixing = time_mixing[sort_idx]
    m_mixing = m_mixing[sort_idx]

    # Drop duplicate times, keeping the last
    if len(time_mixing) > 0:
        diff_mask = np.append(np.diff(time_mixing) != 0, True)
        time_mixing = time_mixing[diff_mask]
        m_mixing = m_mixing[diff_mask]

    min_points = 4 if free_plateau else 3
    if len(time_mixing) < min_points:
        raise ValueError(f"At least {min_points} valid data points required to fit the Lacey curve")

    m0 = m_mixing[0]
    t0 = time_mixing[0]

    m_plateau = None

    if len(m_mixing) >= window_size:
        plateau_window = m_mixing[-window_size:]
        m_plateau = float(np.mean(plateau_window))

        # Warn if tail values are unusually spread around the plateau mean.
        plateau_std = float(np.std(plateau_window))
        if plateau_std > 0:
            lower_bound = m_plateau - (1.0 * plateau_std)
            upper_bound = m_plateau + (1.0 * plateau_std)
            outside_mask = (plateau_window < lower_bound) | (plateau_window > upper_bound)
            n_outside = int(np.sum(outside_mask))
            if n_outside > 0:
                warnings.warn(
                    f"{n_outside} value(s) in the final {window_size} LMI samples fall"
                    f" outside mean +/- 1 std (range: [{lower_bound:.4f}, {upper_bound:.4f}]).",
                    UserWarning,
                )
    
    if m_plateau is None:
        warnings.warn(
            f"Fewer than {window_size} values available. "
            "Using the mean of all available LMI values for m_plateau.",
            UserWarning,
        )
        m_plateau = float(np.mean(m_mixing))

    if free_plateau:
        def fit_func(time, k, tau, m_plat):
            return lacey_mixing_curve_extended(time, k, tau, m0, m_plat)
        
        p0 = (k0, t0, m_plateau)
        bounds = ([0, t0, 1e-10], [np.inf, np.inf, 1.0])
    else:
        fit_func = partial(lacey_mixing_curve_extended, m0=m0, m_plateau=m_plateau)
        p0 = (k0, t0)
        bounds = ([0, t0], [np.inf, np.inf])

    popt, pcov = curve_fit(fit_func, 
                        time_mixing, 
                        m_mixing,
                        p0=p0, 
                        bounds=bounds,
                        maxfev=10000,
                        )
    
    if free_plateau:
        m_plateau = float(popt[2])

    m_fit = fit_func(time_mixing, *popt)
    r2 = r2_score(m_mixing, m_fit)
   
    t95_value = 0.95 * (m_plateau)
    t95_data = None
    above95 = np.where(m_mixing >= t95_value)[0]
    if above95.size > 0:
        t95_data = float(time_mixing[above95[0]])
    t95_fit = popt[1] - (1 / popt[0]) * np.log((0.05 * m_plateau) / (m_plateau - m0))

    return popt, pcov, time_mixing, m_mixing, m_fit, m_plateau, r2, t95_data, t95_fit


# Old functions without m_plateau parameter.
# Kept for analysing any data that does not fully mix within same time frame as bulk of data.


def old_lacey_mixing_curve_extended(time, k, tau, m0):
    """Curve for the Lacey mixing model.

    The Lacey mixing model is a exponential model that describes the 
    mixing of a binary system. The model is classically defined as [1] 
    [2]:

    .. math:: 
    
        M(t) = 1 - (1 - M_0) e^{-kt}

    where :math:`M(t)` is the Lacey mixing index at time :math:`t`, 
    :math:`M_0` is the initial Lacey Mixing index, :math:`k` is the 
    mixing rate constant, and :math:`t` is the time.

    An extension of the Lacey mixing model to include a time constant
    :math:`\\tau` to allow for delayed onset of exponential mixing 
    mixing behavior was proposed by Ratnayake et al. [2]:

    .. math::

        M(t) = max(1 - (1 - M_0) e^{-k(t - \\tau)})

    This mixing model is used to fit the mixing data as with no delayed
    onset of exponential mixing Ratnayake et al.'s model reduces to the
    classical Lacey mixing model.

    References
    ----------

    [1] Lacey PM. Developments in the theory of particle mixing. 
        Journal of applied chemistry. 1954 May;4(5):257-68.

    [2] Ratnayake P, Chandratilleke R, Bao J, Shen Y. A soft-sensor 
        approach to mixing rate determination in powder mixers. Powder 
        Technology. 2018 Aug 1;336:493-505.

    
    Parameters
    ----------
    time : array-like
        The time data for the lacey mixing index.
    k : float
        The rate constant with units of inverse time.
    tau : float
        The time constant. Allows for delayed onset of exponential 
        mixing.
    m0 : float
        The initial value of the Lacey mixing index.

    Returns
    -------
    array-like
        The predicted values for the Lacey mixing curve.

    Raises
    ------
    ValueError
        If time is not an array-like object.
    ValueError
        If k is not an integer or float.
    ValueError
        If tau is not an integer or float.
    ValueError
        If m0 is not an integer or float.
    """
    if not isinstance(time, np.ndarray):
        time = np.array(time)
        raise ValueError("time must be an array-like object")
    
    if not isinstance(k, (int, float)):
        raise ValueError("k must be an integer or float")
    
    if not isinstance(tau, (int, float)):
        raise ValueError("tau must be an integer or float")
    
    if not isinstance(m0, (int, float)):
        raise ValueError("m0 must be an integer or float")
    
    return [max((1 - (1 - m0) * np.exp(-k*(t - tau))), m0) for t in time]

def old_lacey_mixing_curve_fit_extended(time, m, t0=0, tend=None):
    """Fit the Lacey mixing curve to the data.

    The Lacey mixing model is a exponential model that describes the 
    mixing of a binary system. The model is classically defined as [1] 
    [2]:

    .. math:: 
    
        M(t) = 1 - (1 - M_0) e^{-kt}

    where :math:`M(t)` is the Lacey mixing index at time :math:`t`, 
    :math:`M_0` is the initial Lacey Mixing index, :math:`k` is the 
    mixing rate constant, and :math:`t` is the time.

    An extension of the Lacey mixing model to include a time constant
    :math:`\\tau` to allow for delayed onset of exponential mixing 
    mixing behavior was proposed by Ratnayake et al. [2]:

    .. math::

        M(t) = max(1 - (1 - M_0) e^{-k(t - \\tau)})

    This mixing model is used to fit the mixing data as with no delayed
    onset of exponential mixing Ratnayake et al.'s model reduces to the
    classical Lacey mixing model.

    ### MEng-RP additions
    ## Mixing Time, :math:`t_{95}`

    The mixing time of a system can be characterised using the time taken to reach 95%
    of the maximum degree of mixing, :math:`M_{95}`, and is denoted as :math:`t_{95}`. 
    :math:`M_{95}` can be calculated using the parameters from the fit as:

    .. math::
        M_{95} = 0.95

    The :math:`t_{95}` of the input data is the time at which the first value of the 
    data set reaches or exceeds :math:`M_{95}`, if this is not reached then None is returned. 
    The :math:`t_{95}` of the fitted curve is calculated by rearranging the model of t, for a
    value of M_{95}:

    .. math::
        t_{95} = \\tau - \\frac{1}{k} \\ln \\left( \\frac{0.05}{1 - M_0} \\right)

    References
    ----------

    [1] Lacey PM. Developments in the theory of particle mixing. 
        Journal of applied chemistry. 1954 May;4(5):257-68.

    [2] Ratnayake P, Chandratilleke R, Bao J, Shen Y. A soft-sensor 
        approach to mixing rate determination in powder mixers. Powder 
        Technology. 2018 Aug 1;336:493-505.

    
    Parameters
    ----------
    time : array-like
        The time data for the lacey mixing index.
    m : array-like
        The lacey mixing index data.
    t0 : int or float
        The time at which mixing begins by default 0.
    tend : int or float, optional
        The time at which mixing ends, by default None. If None, then
        all the time data from the start time will be used.

    Returns
    -------
    popt : array-like
        The optimal values for the parameters k and tau.
    pcov : 2D array
        The estimated covariance of popt as returned by curve_fit.
    time_mixing : array-like
        The time data for the mixing period used in the fit.
    m_mixing : array-like
        The lacey mixing index data for the mixing period used in the
        fit.
    m_fit : array-like
        The predicted lacey mixing index values for the mixing period 
        calculated using the optimal parameters.
    r2 : float
        The coefficient of determination (R^2) for the fit.
    t95_data : float or None
        The first time in the data at which LMI reaches or exceeds 0.95. None if not reached.
    t95_fit : float
        The time at which the fitted Lacey mixing index reaches/will reach 0.95.

    Raises
    ------
    ValueError
        If time is not an array-like object.
    ValueError
        If m is not an array-like object.
    ValueError
        If t0 is not an integer or float.
    ValueError
        If tend is not an integer or float.
    ValueError
        If time and m are not the same length.
    """
    if not isinstance(time, np.ndarray):
        time = np.array(time)
        raise ValueError("time must be an array-like object")
    
    if not isinstance(m, np.ndarray):
        m = np.array(m)
        raise ValueError("m must be an array-like object")
    
    if not isinstance(t0, (int, float)):
        raise ValueError("t0 must be an integer or float")
    
    if tend is None:
        tend = time[-1]
    elif not isinstance(tend, (int, float)):
        raise ValueError("tend must be an integer or float")
    
    if len(time) != len(m):
        raise ValueError("time and m must be the same length")
    
    mixing_indices = (time >= t0) & (time <= tend)
    time_mixing = time[mixing_indices]
    m_mixing = m[mixing_indices]

    m0 = m_mixing[0]
    t0 = time_mixing[0]

    partial_lacey_mixing_curve = partial(old_lacey_mixing_curve_extended, m0=m0)
    popt, pcov = curve_fit(partial_lacey_mixing_curve, 
                        time_mixing, 
                        m_mixing,
                        p0=(0, t0), 
                        bounds=([0, t0], [np.inf, np.inf]),
                        maxfev=10000,
                        )
    
    m_fit = partial_lacey_mixing_curve(time_mixing, *popt)
    r2 = r2_score(m_mixing, m_fit)

    t95_value = 0.95 
    t95_data = None
    above95 = np.where(m_mixing >= t95_value)[0]
    if above95.size > 0:
        t95_data = float(time_mixing[above95[0]])
    t95_fit = popt[1] - (1 / popt[0]) * np.log((0.05) / (1 - m0))

    return popt, pcov, time_mixing, m_mixing, m_fit, r2, t95_data, t95_fit


def lacey_mixing_curve(time, k, m0, m_plateau=1.0):
    """
    ## Original DEMToolbox function description from Jack Grogan

    ### Curve for the Lacey mixing model.

    The Lacey mixing model is a exponential model that describes the 
    mixing of a binary system. The model is classically defined as [1] 
    [2]:

    .. math:: 
    
        M(t) = 1 - (1 - M_0) e^{-kt}

    where :math:`M(t)` is the Lacey mixing index at time :math:`t`, 
    :math:`M_0` is the initial Lacey Mixing index, :math:`k` is the 
    mixing rate constant, and :math:`t` is the time.


    ### References

    [1] Lacey PM. Developments in the theory of particle mixing. 
        Journal of applied chemistry. 1954 May;4(5):257-68.


    ## MEng-RP additions

    In my MEng-RP a singular value of resolution was used to sample multiple 
    fill heights, as such  the bin volumes vary across runs. This may mean that 
    :math:`M_0 != 0` and varies across runs, and also may mean that 
    :math:`M_{final} != 1`. The above proposed model already uses :math:`M_0` 
    as a parameter, and as such no changes are made in that respect. However, the 
    final value of the Lacey mixing index is taken to be 1 in these models. For 
    the purpose of time-series analysis in my MEng-RP, the term :math:`M_{plateau}` 
    is used to represent the value of Lacey mixing index when the system cannot 
    become more mixed, where :math:`0<M_{plateau}<=1.` As such the Lacey mixing 
    model is expressed as:

    .. math::
    
        M(t) = M_{plateau} - (M_{plateau} - M_0) e^{-kt}.

    
    Parameters
    ----------
    time : array-like
        The time data for the lacey mixing index.
    k : float
        The rate constant with units of inverse time.
    tau : float
        The time constant. Allows for delayed onset of exponential 
        mixing.
    m0 : float
        The initial value of the Lacey mixing index.
    m_plateau : float, optional
        The plateau value of the Lacey mixing index, by default 1.0.

    Returns
    -------
    array-like
        The predicted values for the Lacey mixing curve.

    Raises
    ------
    ValueError
        If time is not an array-like object.
    ValueError
        If k is not an integer or float.
    ValueError
        If m0 is not an integer or float.
    """
    if not isinstance(time, np.ndarray):
        time = np.array(time)
        raise ValueError("time must be an array-like object")
    
    if not isinstance(k, (int, float)):
        raise ValueError("k must be an integer or float")

    if not isinstance(m0, (int, float)):
        raise ValueError("m0 must be an integer or float")
    
    if not isinstance(m_plateau, (int, float)) or not (0 < m_plateau <= 1):
        raise ValueError("m_plateau must be an integer or float between 0 (exclusive) "
        "and 1 (inclusive)")
    
    return [max((m_plateau - (m_plateau - m0) * np.exp(-k*(t))), m0) for t in time]

def lacey_mixing_curve_fit(time, m, window_size=20, t0=0, tend=None, k0=0.1, free_plateau=False):
    """
    ## Original DEMToolbox function description from Jack Grogan

    ### Fit the Lacey mixing curve to the data.

    The Lacey mixing model is a exponential model that describes the 
    mixing of a binary system. The model is classically defined as [1] 
    [2]:

    .. math:: 
    
        M(t) = 1 - (1 - M_0) e^{-kt}

    where :math:`M(t)` is the Lacey mixing index at time :math:`t`, 
    :math:`M_0` is the initial Lacey Mixing index, :math:`k` is the 
    mixing rate constant, and :math:`t` is the time.

    ### References

    [1] Lacey PM. Developments in the theory of particle mixing. 
        Journal of applied chemistry. 1954 May;4(5):257-68.

    ## MEng-RP additions

    In my MEng-RP a singular value of resolution was used to sample 
    multiple fill heights, as such  the bin volumes vary across runs. 
    This may mean that :math:`M_0 != 0` and varies across runs, and also 
    may mean that :math:`M_{final} != 1`. The above proposed models already 
    use :math:`M_0` as a parameter, and as such no changes are made in that 
    respect. However, the final value of the Lacey mixing index is taken to 
    be 1 in these models. For the purpose of time-series analysis in my MEng-RP, 
    the term :math:`M_{plateau}` is used to represent the value of Lacey 
    mixing index when the system cannot become more mixed, where 
    :math:`0<M_{plateau}<=1.` As such the Lacey mixing model is expressed as:

    .. math::
    
        M(t) = M_{plateau} - (M_{plateau} - M_0) e^{-kt},

    Due to the nature of sampling, :math:`M` continues to vary slightly, even after 
    a "stable" degree of mixing has been reached. :math:`M_{plateau}` was calculated 
    as the mean of the final `window_size` values of LMI in the timeseries.
    In my MEng-RP simulations, this was 20, as this was the point when vibration
    stopped.

    ### Mixing Time, :math:`t_{95}`

    The mixing time of a system can be characterised using the time taken to reach 95%
    of the maximum degree of mixing, :math:`M_{95}`, and is denoted as :math:`t_{95}`. 
    :math:`M_{95}` can be calculated using the parameters from the fit as:

    .. math::
        M_{95} = 0.95 * M_{plateau}

    The :math:`t_{95}` of the input data is the time at which the first value of the 
    data set reaches or exceeds :math:`M_{95}`, if this is not reached then None is returned. 
    The :math:`t_{95}` of the fitted curve is calculated by rearranging the model of t, for a
    value of M_{95}:

    .. math::
        t_{95} = - \\frac{1}{k} \\ln \\left( \\frac{0.05 * M_{plateau}}{M_{plateau} - M_0} \\right)


    Parameters
    ----------
    time : array-like
        The time data for the lacey mixing index.
    m : array-like
        The lacey mixing index data.
    window_size : int, optional
        The number of consecutive values to consider when determining if the 
        system has reached the plateau value of the Lacey mixing index, by default 20.
    t0 : int or float
        The time at which mixing begins by default 0.
    tend : int or float, optional
        The time at which mixing ends, by default None. If None, then
        all the time data from the start time will be used.
    k0 : float, optional
        The initial guess for the rate constant k, by default 0.1.
    free_plateau : bool, optional
        If True, m_plateau is treated as a free parameter during fitting, using
        the empirical value as an initial guess. Defaults to False.

    Returns
    -------
    popt : array-like
        The optimal values for the parameters k (and m_plateau if
        free_plateau=True).
    pcov : 2D array or array-like
        The estimated covariance of popt as returned by curve_fit, if free_plateau=True.
        The estimated variance of k if free_plateau=False.
    time_mixing : array-like
        The time data for the mixing period used in the fit.
    m_mixing : array-like
        The lacey mixing index data for the mixing period used in the fit.
    m_fit : array-like
        The predicted lacey mixing index values for the mixing period calculated 
        using the optimal parameters.
    m_plateau : float
        The calculated plateau value of the Lacey mixing index (empirical
        average if free_plateau=False, optimised value if free_plateau=True).
    r2 : float
        The coefficient of determination of the fit.
    t95_data : float or None
        The first time in the actual data where the Lacey mixing index reaches 
        or exceeds 95% of m_plateau. None if not reached.
    t95_fit : float or None
        The interpolated time in the fitted curve where the Lacey mixing index reaches 
        95% of m_plateau.

    Raises
    ------
    ValueError
        If time is not an array-like object.
    ValueError
        If m is not an array-like object.
    ValueError
        If window_size is not an integer > 0.
    ValueError
        If std_threshold is not a float or integer > 0.
    ValueError
        If t0 is not an integer or float.
    ValueError
        If tend is not an integer or float.
    ValueError
        If time and m are not the same length.
    ValueError
        If k0 is not an integer or float.
    ValueError
        If fewer than 3 valid data points remain after preprocessing.

    Warnings
    --------
    UserWarning
        If fewer data points than `window_size`:
        Mean of all available LMI values used for m_plateau.
    UserWarning
        If LMI value in final `window_size` sample falls outside ± 1 standard deviation of mean:
        User message raised.
    """
    if not isinstance(time, np.ndarray):
        time = np.array(time)
        raise ValueError("time must be an array-like object")
    
    if not isinstance(m, np.ndarray):
        m = np.array(m)
        raise ValueError("m must be an array-like object")
    
    if not isinstance(window_size, int) or window_size <= 0:
        raise ValueError("window_size must be an integer > 0")

    if not isinstance(t0, (int, float)):
        raise ValueError("t0 must be an integer or float")
    
    if tend is None:
        tend = time[-1]
    elif not isinstance(tend, (int, float)):
        raise ValueError("tend must be an integer or float")
    
    if len(time) != len(m):
        raise ValueError("time and m must be the same length")
    
    if not isinstance(k0, (int, float)):
        raise ValueError("k0 must be an integer or float")
    
    mixing_indices = (time >= t0) & (time <= tend)
    time_mixing = time[mixing_indices]
    m_mixing = m[mixing_indices]

    # Preprocessing: drop NaNs
    valid_mask = ~np.isnan(m_mixing) & ~np.isnan(time_mixing)
    time_mixing = time_mixing[valid_mask]
    m_mixing = m_mixing[valid_mask]

    # Sort by time
    sort_idx = np.argsort(time_mixing)
    time_mixing = time_mixing[sort_idx]
    m_mixing = m_mixing[sort_idx]

    # Drop duplicate times, keeping the last
    if len(time_mixing) > 0:
        diff_mask = np.append(np.diff(time_mixing) != 0, True)
        time_mixing = time_mixing[diff_mask]
        m_mixing = m_mixing[diff_mask]

    min_points = 4 if free_plateau else 3
    if len(time_mixing) < min_points:
        raise ValueError(f"At least {min_points} valid data points required to fit the Lacey curve")

    m0 = m_mixing[0]
    t0 = time_mixing[0]

    m_plateau = None

    if len(m_mixing) >= window_size:
        plateau_window = m_mixing[-window_size:]
        m_plateau = float(np.mean(plateau_window))

        # Warn if tail values are unusually spread around the plateau mean.
        plateau_std = float(np.std(plateau_window))
        if plateau_std > 0:
            lower_bound = m_plateau - (1.0 * plateau_std)
            upper_bound = m_plateau + (1.0 * plateau_std)
            outside_mask = (plateau_window < lower_bound) | (plateau_window > upper_bound)
            n_outside = int(np.sum(outside_mask))
            if n_outside > 0:
                warnings.warn(
                    f"{n_outside} value(s) in the final {window_size} LMI samples fall"
                    f" outside mean +/- 1 std (range: [{lower_bound:.4f}, {upper_bound:.4f}]).",
                    UserWarning,
                )
    
    if m_plateau is None:
        warnings.warn(
            f"Fewer than {window_size} values available. "
            "Using the mean of all available LMI values for m_plateau.",
            UserWarning,
        )
        m_plateau = float(np.mean(m_mixing))

    if free_plateau:
        def fit_func(time, k, m_plat):
            return lacey_mixing_curve(time, k, m0, m_plat)
        
        p0 = (k0, t0, m_plateau)
        bounds = ([0, t0, 1e-10], [np.inf, np.inf, 1.0])
    else:
        fit_func = partial(lacey_mixing_curve, m0=m0, m_plateau=m_plateau)
        p0 = (k0, t0)
        bounds = ([0, t0], [np.inf, np.inf])

    popt, pcov = curve_fit(fit_func, 
                        time_mixing, 
                        m_mixing,
                        p0=p0, 
                        bounds=bounds,
                        maxfev=10000,
                        )
    
    if free_plateau:
        m_plateau = float(popt[1])

    m_fit = fit_func(time_mixing, *popt)
    r2 = r2_score(m_mixing, m_fit)
   
    t95_value = 0.95 * (m_plateau)
    t95_data = None
    above95 = np.where(m_mixing >= t95_value)[0]
    if above95.size > 0:
        t95_data = float(time_mixing[above95[0]])
    t95_fit = - (1 / popt[0]) * np.log((0.05 * m_plateau) / (m_plateau - m0))

    return popt, pcov, time_mixing, m_mixing, m_fit, m_plateau, r2, t95_data, t95_fit


# Old functions without m_plateau parameter.
# Kept for analysing any data that does not fully mix within same time frame as bulk of data.


def old_lacey_mixing_curve(time, k, m0):
    """Curve for the Lacey mixing model.

    The Lacey mixing model is a exponential model that describes the 
    mixing of a binary system. The model is classically defined as [1] 
    [2]:

    .. math:: 
    
        M(t) = 1 - (1 - M_0) e^{-kt}

    where :math:`M(t)` is the Lacey mixing index at time :math:`t`, 
    :math:`M_0` is the initial Lacey Mixing index, :math:`k` is the 
    mixing rate constant, and :math:`t` is the time.

    References
    ----------

    [1] Lacey PM. Developments in the theory of particle mixing. 
        Journal of applied chemistry. 1954 May;4(5):257-68.

    Parameters
    ----------
    time : array-like
        The time data for the lacey mixing index.
    k : float
        The rate constant with units of inverse time.
    m0 : float
        The initial value of the Lacey mixing index.

    Returns
    -------
    array-like
        The predicted values for the Lacey mixing curve.

    Raises
    ------
    ValueError
        If time is not an array-like object.
    ValueError
        If k is not an integer or float.
    ValueError
        If m0 is not an integer or float.
    """
    if not isinstance(time, np.ndarray):
        time = np.array(time)
        raise ValueError("time must be an array-like object")
    
    if not isinstance(k, (int, float)):
        raise ValueError("k must be an integer or float")

    if not isinstance(m0, (int, float)):
        raise ValueError("m0 must be an integer or float")
    
    return [max((1 - (1 - m0) * np.exp(-k*(t))), m0) for t in time]

def old_lacey_mixing_curve_fit(time, m, t0=0, tend=None):
    """Fit the Lacey mixing curve to the data.

    The Lacey mixing model is a exponential model that describes the 
    mixing of a binary system. The model is classically defined as [1] 
    [2]:

    .. math:: 
    
        M(t) = 1 - (1 - M_0) e^{-kt}

    where :math:`M(t)` is the Lacey mixing index at time :math:`t`, 
    :math:`M_0` is the initial Lacey Mixing index, :math:`k` is the 
    mixing rate constant, and :math:`t` is the time.

    ### MEng-RP additions
    ## Mixing Time, :math:`t_{95}`

    The mixing time of a system can be characterised using the time taken to reach 95%
    of the maximum degree of mixing, :math:`M_{95}`, and is denoted as :math:`t_{95}`. 
    :math:`M_{95}` can be calculated using the parameters from the fit as:

    .. math::
        M_{95} = 0.95

    The :math:`t_{95}` of the input data is the time at which the first value of the 
    data set reaches or exceeds :math:`M_{95}`, if this is not reached then None is returned. 
    The :math:`t_{95}` of the fitted curve is calculated by rearranging the model of t, for a
    value of M_{95}:

    .. math::
        t_{95} = - \\frac{1}{k} \\ln \\left( \\frac{0.05}{1 - M_0} \\right)

    References
    ----------

    [1] Lacey PM. Developments in the theory of particle mixing. 
        Journal of applied chemistry. 1954 May;4(5):257-68.

    Parameters
    ----------
    time : array-like
        The time data for the lacey mixing index.
    m : array-like
        The lacey mixing index data.
    t0 : int or float
        The time at which mixing begins by default 0.
    tend : int or float, optional
        The time at which mixing ends, by default None. If None, then
        all the time data from the start time will be used.

    Returns
    -------
    k : array-like
        The optimal value for the parameter k .
    pvar : 2D array
        The estimated variance of k as returned by curve_fit.
    time_mixing : array-like
        The time data for the mixing period used in the fit.
    m_mixing : array-like
        The lacey mixing index data for the mixing period used in the
        fit.
    m_fit : array-like
        The predicted lacey mixing index values for the mixing period 
        calculated using the optimal parameters.
    r2 : float
        The coefficient of determination (R^2) for the fit.
    t95_data : float or None
        The first time in the data at which LMI reaches or exceeds 0.95. None if not reached.
    t95_fit : float
        The time at which the fitted Lacey mixing index reaches/will reach 0.95.

    Raises
    ------
    ValueError
        If time is not an array-like object.
    ValueError
        If m is not an array-like object.
    ValueError
        If t0 is not an integer or float.
    ValueError
        If tend is not an integer or float.
    ValueError
        If time and m are not the same length.
    """
    if not isinstance(time, np.ndarray):
        time = np.array(time)
        raise ValueError("time must be an array-like object")
    
    if not isinstance(m, np.ndarray):
        m = np.array(m)
        raise ValueError("m must be an array-like object")
    
    if not isinstance(t0, (int, float)):
        raise ValueError("t0 must be an integer or float")
    
    if tend is None:
        tend = time[-1]
    elif not isinstance(tend, (int, float)):
        raise ValueError("tend must be an integer or float")
    
    if len(time) != len(m):
        raise ValueError("time and m must be the same length")
    
    mixing_indices = (time >= t0) & (time <= tend)
    time_mixing = time[mixing_indices]
    m_mixing = m[mixing_indices]

    m0 = m_mixing[0]
    t0 = time_mixing[0]

    partial_lacey_mixing_curve = partial(old_lacey_mixing_curve, m0=m0)
    k, pvar = curve_fit(partial_lacey_mixing_curve, 
                        time_mixing, 
                        m_mixing,
                        p0=(0, t0), 
                        bounds=([0, t0], [np.inf, np.inf]),
                        maxfev=10000,
                        )
    
    m_fit = partial_lacey_mixing_curve(time_mixing, *k)
    r2 = r2_score(m_mixing, m_fit)

    t95_value = 0.95 
    t95_data = None
    above95 = np.where(m_mixing >= t95_value)[0]
    if above95.size > 0:
        t95_data = float(time_mixing[above95[0]])
    t95_fit = - (1 / k) * np.log((0.05) / (1 - m0))

    return k, pvar, time_mixing, m_mixing, m_fit, r2, t95_data, t95_fit