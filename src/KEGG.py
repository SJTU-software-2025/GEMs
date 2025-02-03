import requests
import re
def get_chemical_formula(chemical_name):
    url = f'https://rest.kegg.jp/get/{chemical_name}'
    response = requests.get(url)
    
    if response.status_code == 200:
        data = response.text
        # 使用正则表达式提取化学式
        match = re.search(r'FORMULA\s+([A-Za-z0-9]+)', data)
        if match:
            formula = match.group(1)
        else:
            print(f"{chemical_name}:未找到化学式。")
            return None
        return formula
    elif response.status_code == 404:
        print(f"404 not found")
        return None
    else:
        print(f"{chemical_name}:未找到化学式。")
        return None

def main():
    # 例子: 获取ATP的分子式
    chemical_name = "C00002"  # KEGG中的ATP ID
    formula = get_chemical_formula(chemical_name)
    print(f"The chemical formula of ATP is: {formula}")

if __name__ == "__main__":
    main()
