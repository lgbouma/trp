import numpy as np
import matplotlib.pyplot as plt
from trp.lcprocessing import calculate_lsp

def test_lsp_calc():
    # Set the parameters for the simulated timeseries
    period = 2.0  # Period in days
    amplitude = 1.0
    A_2 = 0.01
    phase = 0.0
    mean = 0.0
    std_dev = 0.2
    num_points = int(27*48*(1800/200))
    time_range = 27

    # Generate the time values
    times = np.linspace(0, time_range, num_points)

    # Generate the almost-sinusoid signal
    signal = (
        amplitude * np.sin(2 * np.pi * (times / period) + phase) + mean
        + A_2 * np.sin(2 * np.pi * (times / (2*period)) + phase)
    )

    # Add Gaussian noise to the signal
    noise = np.random.normal(0, std_dev, num_points)
    fluxes = signal + noise

    # Calculate the Lomb Scargle periodogram
    d = calculate_lsp(times, fluxes, starid='whateversecret', outdir='./',
                      cachedict={'hi':42})

    # Get the peak period from the d
    peak_period = d['period']

    # Set the tolerance for period comparison
    tolerance = 0.1

    # Assert that the recovered peak period is close to the true period
    assert np.abs(peak_period - period) < tolerance, (
        f"Expected period: {period}, Recovered period: {peak_period}"
    )
    if A_2 <= 0.01:
        assert np.abs(d['reduced_chi2'] - 1) < tolerance, (
            f"Expected reduced_chi2 of 1, got {d['reduced_chi2']}"
        )

    print(d)

    # Print the top N best periods and values
    n_best = len(d['lsp']['nbestperiods'])
    print(f"Top {n_best} best periods and values:")
    for i in range(n_best):
        print(
            f"Period: {d['lsp']['nbestperiods'][i]:.3f}, "
            f"Value: {d['lsp']['nbestlspvals'][i]:.3f}"
        )

    fig, axs = plt.subplots(nrows=2)
    ax = axs[0]
    ax.scatter(times, fluxes, zorder=2, s=0.1)
    ax.plot(times, d['best_fit_sinusoid'], zorder=3, lw=0.5, c='k')
    ax = axs[1]
    ax.plot(d['lsp']['periods'], d['lsp']['power'])
    ax.vlines(d['period'], 0, 1, alpha=1, ls="-", lw=0.5)
    ax.vlines(d['lsp']['nbestperiods'], 0, 1, alpha=0.5, ls=":", lw=0.5)
    ax.set_yscale('log')
    fig.savefig('secret_temp.png')

if __name__ == "__main__":
    test_lsp_calc()
