# healthy diet score by U.S. Dietary Guidelines

Calculate a 0–7 healthy diet score using seven components (vegetables,
fruit, fish, whole grains, refined grains, processed meat, unprocessed
meat) following the U.S. Dietary Guidelines targets. Each component
meeting its intake target contributes 1 point; higher scores indicate
healthier patterns.

## Usage

``` r
leo_diet_us(
  fresh_fruit,
  dried_fruit,
  cooked_veg,
  salad_veg,
  oily_fish,
  non_oily_fish,
  bread_intake,
  bread_type,
  cereal_intake,
  cereal_type,
  processed_meat,
  age_last_meat,
  poultry,
  beef,
  lamb,
  pork,
  ...
)
```

## Arguments

- fresh_fruit:

  Fresh fruit pieces per day.

- dried_fruit:

  Dried fruit pieces per day (5 pieces = 1 serving).

- cooked_veg:

  Cooked vegetables tablespoons per day (3 tbsp = 1 serving).

- salad_veg:

  Salad/raw vegetables tablespoons per day (3 tbsp = 1 serving).

- oily_fish:

  Oily fish frequency code (p1329_i0, 0–5 per Data-Coding 100377).

- non_oily_fish:

  Non-oily fish frequency code (p1339_i0, 0–5 per Data-Coding 100377).

- bread_intake:

  Bread slices per week (p1438_i0); -10 treated as 0/week (less than
  once).

- bread_type:

  Bread type code (p1448_i0): 1=White, 2=Brown, 3=Wholemeal, 4=Other.

- cereal_intake:

  Cereal bowls per week (p1458_i0); -10 treated as 0/week.

- cereal_type:

  Cereal type code (p1468_i0): 1=Bran, 2=Biscuit, 3=Oat, 4=Muesli,
  5=Other.

- processed_meat:

  Processed meat frequency code (p1349_i0, 0–5 per Data-Coding 100377).

- age_last_meat:

  Age when last ate meat (p3680_i0); 0 = never ate meat.

- poultry:

  Poultry frequency code (p1359_i0, 0–5 per Data-Coding 100377).

- beef:

  Beef frequency code (p1369_i0, 0–5 per Data-Coding 100377).

- lamb:

  Lamb/mutton frequency code (p1379_i0, 0–5 per Data-Coding 100377).

- pork:

  Pork frequency code (p1389_i0, 0–5 per Data-Coding 100377).

- ...:

  Reserved for future use.

## Value

Integer vector in 0,7; NA if any component cannot be evaluated for a
row.

## Details

Field expectations (UKB instance 0):

- Fruits: p1309_i0 (fresh pieces/day), p1319_i0 (dried pieces/day; 5
  pieces = 1 serving)

- Vegetables: p1289_i0 (cooked tablespoons/day), p1299_i0 (salad
  tablespoons/day; 3 tbsp = 1 serving)

- Fish: p1329_i0 (oily fish frequency code), p1339_i0 (non-oily fish
  frequency code)

- Bread: p1438_i0 (slices/week, -10 = "less than once"), p1448_i0 (type:
  1=White, 2=Brown, 3=Wholemeal, 4=Other)

- Cereal: p1458_i0 (bowls/week, -10 = "less than once"), p1468_i0 (type:
  1=Bran, 2=Biscuit, 3=Oat, 4=Muesli, 5=Other)

- Processed meat: p1349_i0 (frequency code 0–5, per Data-Coding 100377);
  p3680_i0 == 0 implies never ate meat

- Unprocessed meat: p1359_i0 (poultry), p1369_i0 (beef), p1379_i0
  (lamb), p1389_i0 (pork) — all frequency codes 0–5

All frequency fields (100377) encode: 0=Never, 1=\<1/wk, 2=1/wk,
3=2–4/wk, 4=5–6/wk, 5=Daily. To meet component thresholds, frequency
codes must be converted to weekly servings (e.g., code 2 → 1/wk, code 3
→ 3/wk midpoint, code 5 → 7/wk).

Whole grains: bread_type==3 (wholemeal) + cereal_type in (1,3,4)
(bran/oat/muesli) Refined grains: bread_type in (1,2,4)
(white/brown/other) + cereal_type in (2,5) (biscuit/other)
