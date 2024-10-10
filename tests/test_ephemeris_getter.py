from trp.getters import get_ephemeris_given_ticid

def main():
    # wasp-4
    tic_id = '402026209'
    ephem_df = get_ephemeris_given_ticid(tic_id)
    assert abs(ephem_df['Period'].iloc[0] - 1.338) < 1e-3

if __name__ == "__main__":
    main()
