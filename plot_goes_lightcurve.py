import argparse
import pandas as pd
import matplotlib.pyplot as plt

def plot_goes_lightcurve(csv_file, save_plot=False, output_file="goes_lightcurve.png"):
    # Load CSV
    df = pd.read_csv(csv_file, parse_dates=['Time'], index_col='Time')

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(df.index, df['xrsa'], label='0.5–4 Å (short)', color='red')
    plt.plot(df.index, df['xrsb'], label='1–8 Å (long)', color='blue')

    plt.yscale('log')
    plt.xlabel('Time (UTC)')
    plt.ylabel('Flux (W/m²)')
    plt.title('GOES Soft X-ray Light Curve')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    if save_plot:
        plt.savefig(output_file,dpi=300)
        print(f"Plot saved to {output_file}")
    else:
        plt.show()

def main():
    parser = argparse.ArgumentParser(description="Plot GOES X-ray light curve from CSV.")
    parser.add_argument("csv_file", help="Path to the GOES CSV file")
    parser.add_argument("--save", action="store_true", help="Save plot to file instead of showing it")
    parser.add_argument("--output", default="goes_lightcurve.png", help="Output filename for the plot")

    args = parser.parse_args()

    plot_goes_lightcurve(args.csv_file, save_plot=args.save, output_file=args.output)

if __name__ == "__main__":
    main()
