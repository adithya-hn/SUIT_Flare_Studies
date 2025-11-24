import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
# Load the combined report file
report_file = "combined_flare_report.csv"  # Update this with the actual filename

df = pd.read_csv(report_file)

# Extract unique filters and parameters to plot
filters = df['Filter'].unique()
parameters = [
              'Flare_Contrast',  'Enhancement_QS']

# Set Seaborn style
sns.set_theme(style="whitegrid", palette="Set2")

# Create plots for each parameter
for param in parameters:
    plt.figure(figsize=(10, 6))
    ax = sns.barplot(data=df, x="Flare_ID", y=param, hue="Filter")
    
    # Add value labels on top of bars
    for p in ax.patches:
        ax.annotate(f'{p.get_height():.2f}', 
                    (p.get_x() + p.get_width() / 2., p.get_height()),
                    ha='center', va='bottom', fontsize=10, color='black', fontweight='bold')
    
    plt.xlabel("Flare Case")
    plt.ylabel(param.replace("_", " "))
    plt.title(f"Variation of {param.replace('_', ' ')} across Flare Cases")
    plt.xticks(rotation=45)
    plt.legend(title="Filter")
    plt.tight_layout()
    
    # Save the plot
    plot_filename = f"{param}_variation.png"
    plt.savefig(plot_filename,dpi=300)
    plt.close()
    print(f"Plot saved as {plot_filename}")
