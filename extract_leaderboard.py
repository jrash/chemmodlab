"""
Script to extract leaderboard table from PolarisHub competition page.
Multiple methods provided - try them in order.
"""

import pandas as pd
import json

# Method 1: Using Selenium (for JavaScript-rendered content)
def extract_with_selenium(url):
    """
    Best method for dynamically loaded content.
    Requires: pip install selenium webdriver-manager
    """
    try:
        from selenium import webdriver
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import WebDriverWait
        from selenium.webdriver.support import expected_conditions as EC
        from selenium.webdriver.chrome.options import Options
        from webdriver_manager.chrome import ChromeDriverManager
        from selenium.webdriver.chrome.service import Service

        # Set up headless Chrome
        chrome_options = Options()
        chrome_options.add_argument("--headless")
        chrome_options.add_argument("--no-sandbox")
        chrome_options.add_argument("--disable-dev-shm-usage")

        driver = webdriver.Chrome(
            service=Service(ChromeDriverManager().install()),
            options=chrome_options
        )

        print("Loading page...")
        driver.get(url)

        # Wait for table to load
        wait = WebDriverWait(driver, 10)
        wait.until(EC.presence_of_element_located((By.TAG_NAME, "table")))

        # Extract table data
        tables = driver.find_elements(By.TAG_NAME, "table")

        if tables:
            # Get the first table (or iterate through all)
            table = tables[0]

            # Extract headers
            headers = []
            header_row = table.find_elements(By.TAG_NAME, "th")
            for th in header_row:
                headers.append(th.text)

            # Extract rows
            rows_data = []
            rows = table.find_elements(By.TAG_NAME, "tr")
            for row in rows[1:]:  # Skip header row
                cols = row.find_elements(By.TAG_NAME, "td")
                row_data = [col.text for col in cols]
                if row_data:
                    rows_data.append(row_data)

            # Create DataFrame
            df = pd.DataFrame(rows_data, columns=headers)
            driver.quit()
            return df
        else:
            driver.quit()
            print("No table found on page")
            return None

    except ImportError:
        print("Selenium not installed. Run: pip install selenium webdriver-manager")
        return None
    except Exception as e:
        print(f"Error with Selenium: {e}")
        return None


# Method 2: Using requests and BeautifulSoup
def extract_with_bs4(url):
    """
    Works for static HTML tables.
    Requires: pip install requests beautifulsoup4 lxml
    """
    try:
        import requests
        from bs4 import BeautifulSoup

        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
        }

        response = requests.get(url, headers=headers)
        soup = BeautifulSoup(response.content, 'lxml')

        # Try pandas read_html first
        tables = pd.read_html(str(soup))
        if tables:
            return tables[0]

        return None

    except ImportError:
        print("BeautifulSoup not installed. Run: pip install requests beautifulsoup4 lxml")
        return None
    except Exception as e:
        print(f"Error with BeautifulSoup: {e}")
        return None


# Method 3: Manual extraction using browser DevTools
def print_manual_instructions():
    """
    Instructions for manual extraction using browser DevTools
    """
    instructions = """
    MANUAL EXTRACTION METHOD:

    1. Open the page in your browser
    2. Open Developer Tools (F12 or Right-click -> Inspect)
    3. Go to Console tab
    4. Paste this JavaScript code:

    // Extract table data
    let table = document.querySelector('table');
    let headers = Array.from(table.querySelectorAll('th')).map(th => th.innerText);
    let rows = Array.from(table.querySelectorAll('tr')).slice(1).map(row => {
        return Array.from(row.querySelectorAll('td')).map(td => td.innerText);
    });

    let data = rows.map(row => {
        let obj = {};
        headers.forEach((header, i) => {
            obj[header] = row[i];
        });
        return obj;
    });

    console.log(JSON.stringify(data, null, 2));
    copy(JSON.stringify(data, null, 2));

    5. The data will be copied to clipboard as JSON
    6. Save it to a file called 'leaderboard_data.json'
    7. Then run: load_from_json('leaderboard_data.json')
    """
    print(instructions)


def load_from_json(filepath):
    """Load data from manually extracted JSON"""
    try:
        with open(filepath, 'r') as f:
            data = json.load(f)
        df = pd.DataFrame(data)
        return df
    except Exception as e:
        print(f"Error loading JSON: {e}")
        return None


# Method 4: Using pandas directly
def extract_with_pandas(url):
    """
    Simplest method - works if table is in simple HTML
    """
    try:
        tables = pd.read_html(url)
        if tables:
            print(f"Found {len(tables)} table(s)")
            return tables[0]  # Return first table
        return None
    except Exception as e:
        print(f"Error with pandas: {e}")
        return None


def main():
    url = "https://polarishub.io/competitions/asap-discovery/antiviral-drug-discovery-2025"

    print("Attempting to extract leaderboard table...\n")

    # Try methods in order of likelihood to work
    methods = [
        ("Selenium (JavaScript rendering)", extract_with_selenium),
        ("Pandas (simple)", extract_with_pandas),
        ("BeautifulSoup", extract_with_bs4),
    ]

    df = None
    for method_name, method_func in methods:
        print(f"\nTrying: {method_name}")
        df = method_func(url)
        if df is not None:
            print(f"✓ Success with {method_name}!")
            break

    if df is not None:
        # Save to CSV
        output_file = "leaderboard_data.csv"
        df.to_csv(output_file, index=False)
        print(f"\n✓ Data saved to: {output_file}")
        print(f"\nPreview:\n{df.head()}")
        print(f"\nShape: {df.shape}")
        return df
    else:
        print("\n⚠ Automated methods failed. Trying manual method:")
        print_manual_instructions()
        return None


if __name__ == "__main__":
    df = main()
