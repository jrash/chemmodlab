# Extracting the PolarisHub Leaderboard Table

Multiple methods to extract the leaderboard table from:
`https://polarishub.io/competitions/asap-discovery/antiviral-drug-discovery-2025`

## Method 1: Browser Console (EASIEST - RECOMMENDED) ⭐

This is the quickest method that works for any dynamically loaded table:

1. Open the competition page in your browser
2. Press `F12` (or right-click → Inspect) to open Developer Tools
3. Click on the **Console** tab
4. Copy the entire contents of `browser_extract.js` and paste into the console
5. Press Enter

The script will:
- Automatically extract all table data
- Copy JSON to your clipboard
- Display the data in console
- Provide download functions

To download the extracted data:
```javascript
downloadAsJSON()  // Downloads leaderboard.json
downloadAsCSV()   // Downloads leaderboard.csv
```

## Method 2: Python with Selenium (AUTOMATED)

Best for automated/repeated extractions:

### Setup:
```bash
pip install -r requirements_extraction.txt
```

### Run:
```bash
python extract_leaderboard.py
```

This will try multiple extraction methods automatically and save to `leaderboard_data.csv`.

## Method 3: Manual Copy-Paste

If you just need the data once:

1. Open the page in your browser
2. Select the entire table
3. Copy (Ctrl+C / Cmd+C)
4. Paste into Excel, Google Sheets, or a text editor
5. Save as CSV

## Expected Table Structure

The leaderboard table should contain columns like:
- **Name**: Model/submission name
- **Contributors**: Team members
- **CLD**: Chemical space coverage metric
- **Mean Absolute Error**: MAE metric
- **Mean Squared Error**: MSE metric
- **Pearson R**: Pearson correlation coefficient
- **Spearman Rho**: Spearman rank correlation
- **R²**: Coefficient of determination
- **Kendall τ**: Kendall tau correlation
- **References**: Links to code/papers

## Troubleshooting

### "No table found"
- Make sure the page is fully loaded before running the script
- The table might be in a different tab (check Task dropdown for "MLM")
- Try scrolling down to ensure the table is rendered

### Python script fails
- The website might be blocking automated access (403 error)
- Try the browser console method instead
- Or use a headless browser with Selenium

### Data looks incomplete
- The table might be paginated - check for "Previous/Next" buttons
- You may need to click "Next" and run the extraction multiple times
- Or modify the script to handle pagination

## Using the Extracted Data

Once you have the CSV file:

```python
import pandas as pd

# Load the data
df = pd.read_csv('leaderboard_data.csv')

# View summary
print(df.head())
print(df.describe())

# Sort by metric
df_sorted = df.sort_values('Mean Absolute Error')

# Filter top performers
top_5 = df.nsmallest(5, 'Mean Absolute Error')
```

## Notes

- The table appears to be dynamically loaded with JavaScript
- Browser console method is most reliable
- Data structure may change as competition progresses
- Respect the website's terms of service when scraping
