import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Data: median ± std for each method at different i values
data = {
    5: {
        'gamma & gaussian': (39.66, 0.74),
        'skew normal 1S2D': (26.32, 2.84),
        'skew normal 2S3D': (37.62, 3.19),
    },
    7: {
        'gamma & gaussian': (31.84, 0.62),
        'skew normal 1S2D': (27.78, 0.16),
        'skew normal 2S3D': (31.27, 0.13),
    },
    9: {
        'gamma & gaussian': (34.40, 0.40),
        'skew normal 1S2D': (29.30, 0.08),
        'skew normal 2S3D': (31.47, 0.10),
    },
    11: {
        'gamma & gaussian': (36.05, 0.15),
        'skew normal 1S2D': (33.50, 0.41),
        'skew normal 2S3D': (44.39, 0.31),
    },
    14: {
        'gamma & gaussian': (35.01, 0.11),
        'skew normal 1S2D': (31.68, 0.17),
        'skew normal 2S3D': (44.47, 0.16),
    },
    18: {
        'gamma & gaussian': (44.48, 0.21),
        'skew normal 1S2D': (37.42, 1.73),
        'skew normal 2S3D': (31.28, 0.12),
    },
    21: {
        'gamma & gaussian': (33.78, 0.22),
        'skew normal 1S2D': (30.71, 0.88),
        'skew normal 2S3D': (37.98, 5.19),
    },
    25: {
        'gamma & gaussian': (44.79, 0.49),
        'skew normal 1S2D': (36.21, 0.21),
        'skew normal 2S3D': (32.87, 0.31),
    },
    30: {
        'gamma & gaussian': (35.33, 0.22),
        'skew normal 1S2D': (31.15, 0.09),
        'skew normal 2S3D': (31.41, 0.19),
    },
}

# Create output directory if it doesn't exist
output_dir = './bootstrap_boxplots'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Loop through each i value
for i, methods in data.items():
    # Extract the values for the three methods
    gamma_gaussian = np.random.normal(loc=methods['gamma & gaussian'][0], scale=methods['gamma & gaussian'][1], size=100)
    skew_normal_1S2D = np.random.normal(loc=methods['skew normal 1S2D'][0], scale=methods['skew normal 1S2D'][1], size=100)
    skew_normal_2S3D = np.random.normal(loc=methods['skew normal 2S3D'][0], scale=methods['skew normal 2S3D'][1], size=100)

    # Combine the data for the boxplot
    data_to_plot = [gamma_gaussian, skew_normal_1S2D, skew_normal_2S3D]

    # Calculate sigma squared (standard deviation squared)
    sigma_squared = [
        methods['gamma & gaussian'][1]**2,
        methods['skew normal 1S2D'][1]**2,
        methods['skew normal 2S3D'][1]**2
    ]
    
    # Create the boxplot
    plt.figure(figsize=(8, 6))
    sns.boxplot(data=data_to_plot)
    
    # Set the title, labels, and x-ticks
    plt.title(f'Comparison of 3 Methods for i = {i}')
    plt.xlabel('Method')
    plt.ylabel('Value')
    plt.xticks([0, 1, 2], ['Gamma & Gaussian', 'Skew Normal 1S2D', 'Skew Normal 2S3D'])
    
    # Annotate the sigma squared (σ²) values on the plot
    for j, sigma_sq in enumerate(sigma_squared):
        plt.text(j, max(data_to_plot[j]) + 0.5, f'$\sigma^2$ = {sigma_sq:.2f}', 
                 horizontalalignment='center', fontsize=10, color='red')

    # Save the plot as a PNG file
    plt.savefig(os.path.join(output_dir, f'synthetic_{i}.png'))
    
    # Close the plot to free up memory
    plt.close()

print("Boxplots saved in './bootstrap_boxplots'")

