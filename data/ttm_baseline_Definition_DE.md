# Definition des Datensatzes `ttm_baseline.csv`

**Zweck.** `ttm_baseline.csv` ist der pro-Patient-Basisdatensatz für die
Überlebenszeitanalyse der malignen Transformation (MT) bei niedriggradigen
Gliomen (LGG). Eine Zeile pro Patient. Er bildet das intervallzensierte
Ereigniszeit-Outcome ab und enthält die klinischen, molekularen und
bildgebenden Basiskovariaten.

**Umfang.** 155 Patienten · 12 Spalten · 77 MT-Ereignisse · 78 zensiert.

---

## Herkunft der Variablen (wichtig)

Der Datensatz wurde aus **zwei** Quellen zusammengeführt; das erklärt, warum
einige Spalten in den `*.rds`-Dateien fehlen:

| Quelle | Datei | liefert |
|---|---|---|
| Präoperativer Klinik-/Radiomics-Export | `Clinical_preop.rds` | `patient_number`, `age_atDiag45`, `vol_diag_centered`, `contrast_enhancingI` |
| **Klinische Rohdatenbank** | `Data_LGG_Clinical_update.xlsx`, Blatt „Datensaetzte" | `idh`, `p19q`, `CDKAN2A/B`, `surgerytype` (→ `EOR`), Operations-, MRT- und Nachbeobachtungs-Zeitpunkte (→ `L`, `R`, `event`, `case`) |

> **Hinweis zur Rückfrage der Statistikerin:** Die Variablen `EOR`, `idh`,
> `p19q` und `CDKAN2A/B` sind in **keiner** der `*.rds`-Dateien enthalten —
> sie stammen ausschließlich aus `Data_LGG_Clinical_update.xlsx`.
> `contrast_enhancingI` ist hingegen sehr wohl in `Clinical_preop.rds`
> vorhanden. `L` und `R` sind **keine** direkte Kopie der RDS-Spalten
> `time_dia_mali_1OP_low/up`, sondern werden aus den Klinik-Zeitpunkten neu
> konstruiert (siehe unten).

Die Datei `Data_LGG_Clinical_update.xlsx` enthält identifizierende Daten
(Geburtsdatum, OP-Daten, Vitalstatus) und wird daher institutionell unter
EK-2021-1042 verwahrt; sie liegt **nicht** im öffentlichen GitHub-Repository.

---

## Variablendefinitionen

| Spalte | Typ | Definition |
|---|---|---|
| `patient_number` | ganzzahlig | Anonymisierte Patienten-ID (Schlüssel über alle Datensätze). |
| `L` | numerisch | **Untere Grenze des MT-Intervalls**, in Jahren ab Diagnose. Siehe Intervall-Logik unten. |
| `R` | numerisch | **Obere Grenze des MT-Intervalls**, in Jahren ab Diagnose; `Inf` bei rechtszensierten Patienten. |
| `event` | 0 / 1 | 1 = MT beobachtet (R endlich); 0 = zensiert (R = `Inf`). |
| `case` | Case_0 / Case_2 / Case_3 | Fallklassifikation der MT-Beobachtung (s. u.). |
| `EOR` | Faktor | Resektionsausmaß: GTR / STR / Biopsy. Rekodiert aus `surgerytype` (s. u.). |
| `idh` | 0 / 1 | IDH-Mutationsstatus (1 = mutiert). Quelle: Klinik-DB. |
| `p19q` | 0 / 1 | 1p/19q-Kodeletion (1 = vorhanden). Quelle: Klinik-DB. |
| `CDKAN2A/B` | Text | CDKN2A/B-Alteration: „No" / „Yes" / „Yes (heterozygot)". Quelle: Klinik-DB. |
| `age_atDiag45` | numerisch | Alter bei Diagnose, **zentriert bei 45 Jahren** (Rohalter = Wert + 45). |
| `vol_diag_centered` | numerisch | Tumorvolumen bei Diagnose in cm³, **median-zentriert**. |
| `contrast_enhancingI` | „N" / „Y" | Kontrastmittelanreicherung in der Baseline-Bildgebung (Y = anreichernd). |

---

## Intervallzensierungs-Logik (`L`, `R`, `event`, `case`)

Die MT-Ereigniszeit ist **intervallzensiert**: das tatsächliche
Transformationsdatum liegt zwischen zwei aufeinanderfolgenden MRT-Studien. Das
Intervall `[L, R]` wird patientenspezifisch aus den Klinik-Zeitpunkten
konstruiert — **nicht** aus den RDS-Spalten `time_dia_mali_1OP_low/up` kopiert.

| `case` | n | Bedeutung | `L` | `R` | `event` |
|---|---|---|---|---|---|
| **Case_0** | 78 | Keine MT beobachtet → rechtszensiert | letzter MRT-/Nachbeobachtungszeitpunkt | `Inf` | 0 |
| **Case_2** | 27 | MT **vor** der ersten Operation beobachtet | 0 (Diagnose) | erste MRT mit MT-Nachweis (bzw. Transformations-OP) | 1 |
| **Case_3** | 50 | MT **nach** der Operation beobachtet | letzte MRT ohne Malignitätszeichen | erste MRT mit MT-Nachweis | 1 |

Damit gilt durchgängig: `event = 1 ⇔ R` endlich; `event = 0 ⇔ R = Inf`.

---

## Rekodierung von `EOR` aus `surgerytype`

`EOR` ist **keine** eigenständige Spalte einer Quelldatei, sondern wird aus
`surgerytype` der Klinik-DB rekodiert (Kodierung laut Blatt „Kodierung"):

| `surgerytype` | Klartext | `EOR` | n |
|---|---|---|---|
| 1 | Resection | **GTR** (gross total resection) | 84 |
| 3 | Subtotal Resection | **STR** | 47 |
| 2 | Biopsy | **Biopsy** | 24 |

(Validiert: Kreuztabelle `EOR × surgerytype` ist exakt 1:1, ohne Mehrdeutigkeit.)

---

## Validierungsstand

- `idh`, `p19q`, `CDKAN2A/B` aus `ttm_baseline.csv` stimmen mit
  `Data_LGG_Clinical_update.xlsx` für alle 155 Patienten überein.
- `EOR ↔ surgerytype` ist eine exakte 1:1-Abbildung.
- `event ↔ Endlichkeit von R` ist konsistent (77 Ereignisse / 78 zensiert).
- Radiomics (`radiomics_99_zscore.csv`, `radiomics_428_zscore.csv`) sind
  bereits **z-standardisiert** (Mittelwert 0, SD 1) und stimmen bis auf
  Gleitkomma-Rundung (≈ 5·10⁻¹⁵) mit `Radiomics_99.rds` / `Radiomics_428.rds`
  überein.

---

*Erstellt zur Klärung der Variablen-Herkunft für die statistische
Reproduktion. Die L/R-Konstruktion und der Klinik→CSV-Zusammenführungsschritt
sind nicht als Skript im Repository hinterlegt; die zugrunde liegenden
Zeitpunkte liegen in `Data_LGG_Clinical_update.xlsx`.*
