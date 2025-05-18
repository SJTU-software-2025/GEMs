https://www.nature.com/articles/s41467-020-20338-2
(https://www.researchsquare.com/article/rs-1823493/v1)

https://www.nature.com/articles/s41598-023-32982-x

https://scholar.google.com/scholar?hl=zh-CN&as_sdt=0%2C5&q=etfl&btnG=


---

https://pmc.ncbi.nlm.nih.gov/articles/PMC11033187/

总酶量约束公式的推导基于酶动力学参数和细胞蛋白质组分配的限制，其核心思想是确保代谢模型中所有反应所需的酶总量不超过细胞可提供的酶能力。

\[
\sum_{i=1}^{n} \left( v_i \cdot MW_i \cdot \frac{\sigma_i}{k_{\text{cat},i}} \right) \leq P_{\text{total}} \cdot f
\]
其中：

- \(v_i\)：反应 \(i\) 的通量（单位：mmol/gDCW/h）；
- \(MW_i\)：催化反应 \(i\) 的酶的分子量（单位：g/mol）；
- \(\sigma_i\)：酶的饱和系数（无量纲，通常取0.5）；
- \(k_{\text{cat},i}\)：酶的转换数（单位：s⁻¹）；
- \(P_{\text{total}}\)：总细胞蛋白质含量（单位：g/gDCW，如0.56 g/gDCW）；
- \(f\)：酶占总蛋白质的质量分数（如大肠杆菌 \(f = 0.406\)）。

**推导步骤**

**1. 单反应所需酶量**

每个反应 \(i\) 的酶需求由以下因素决定：

- **催化效率**：酶的 \(k_{\text{cat}}\) 表示单个酶分子每秒催化的底物摩尔数。
- **通量需求**：反应通量 \(v_i\) 转换为每秒单位（1 h = 3600 s）：
  \[
  v_i \, (\text{mmol/gDCW/h}) = \frac{v_i}{3600} \, (\text{mmol/gDCW/s}).
  \]
- **酶浓度**：为维持通量 \(v_i\)，所需酶浓度为：
  \[
  E_i = \frac{v_i}{k_{\text{cat},i} \cdot 3600} \, (\text{mmol/gDCW}).
  \]

**2. 酶质量计算**

将酶浓度转换为质量（考虑分子量 \(MW_i\)）：
\[
E_i \, (\text{g/gDCW}) = \frac{v_i}{k_{\text{cat},i} \cdot 3600} \cdot MW_i \cdot 10^{-3}.
\]
（注：\(10^{-3}\) 转换 mmol 为 mol。）

**3. 饱和系数修正**

酶在实际条件下可能未完全饱和（如底物不足或存在抑制），引入饱和系数 \(\sigma_i\)（通常取0.5）：
\[
E_i^{\text{实际}} = E_i \cdot \frac{1}{\sigma_i}.
\]

**4. 总酶量约束**

所有反应的总酶需求需满足：
\[
\sum_{i=1}^{n} E_i^{\text{实际}} \leq P_{\text{total}} \cdot f.
\]
代入后得到最终公式：
\[
\sum_{i=1}^{n} \left( \frac{v_i \cdot MW_i}{k_{\text{cat},i} \cdot \sigma_i \cdot 3600} \cdot 10^{-3} \right) \leq P_{\text{total}} \cdot f.
\]
**简化形式**（省略单位转换常数）：
\[
\sum_{i=1}^{n} \left( v_i \cdot MW_i \cdot \frac{\sigma_i}{k_{\text{cat},i}} \right) \leq P_{\text{total}} \cdot f \cdot 3600 \cdot 10^{3}.
\]
实际应用中，单位转换常数通常隐含在参数定义中，因此直接使用原始公式。

**关键参数来源**

1. **饱和系数（\(\sigma_i\)）**：  
   通常假设为0.5，基于实验观测（如酶平均活性为最大活性的50%）。
2. **总蛋白质含量（\(P_{\text{total}}\)）**：  
   微生物细胞中平均为0.56 g/gDCW。
3. **酶质量分数（\(f\)）**：  
   不同物种差异显著（如大肠杆菌 \(f=0.406\)，酵母 \(f=0.446\)）。

