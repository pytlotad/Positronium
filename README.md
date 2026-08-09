# Klasyczne pozytonium — symulacja elektrodynamiczna

Program w CERN ROOT przedstawia klasyczny model elektrodynamiczny układu
elektron–pozyton:

- przyciąganie cząstek jest liczone w kontrolowanym, niskoprędkościowym
  przybliżeniu Coulomba (początkowo `v/c ≈ 0,005`);
- pole magnetyczne poruszających się ładunków jest uwzględniane w przybliżeniu
  niskich prędkości wraz z siłą Lorentza `q v × B`;
- ruch jest całkowany dla relatywistycznego pędu `p = γmv`;
- reakcja promieniowania wynika z nierelatywistycznej redukcji rzędu
  Landaua–Lifshitza `F_RR = τ dF_ext/dt`; pochodna całej siły zewnętrznej jest
  liczona symetrycznie wzdłuż chwilowej trajektorii;
- moc promieniowania jest niezależnie liczona wzorem Liénarda i służy do
  kontroli wyniku, ale nie jest używana do korygowania trajektorii.

Każda cząstka ma też klasyczny moment dipolowy o losowym kierunku początkowym.
Losowane są również radialna i styczna prędkość względna. Wybrane w menu
zjawisko określa kontrolowany zakres losowania prowadzący do bezpośredniego
zderzenia, rozproszenia albo związanego stanu para-/orto-podobnego. Momenty
dipolowe zachowują w trakcie pojedynczej symulacji wylosowaną orientację.
Siła dipol–dipol i jej energia potencjalna są włączone do ruchu orbitalnego,
choć przy promieniu Bohra pozostają wiele rzędów wielkości słabsze od
przyciągania Coulomba. Spin pozytonium jest w istocie zjawiskiem kwantowym,
więc dipole stanowią tutaj jedynie klasyczny model fenomenologiczny.
Punktowe oddziaływania magnetyczne są gładko wygaszane poniżej `0,7 a₀`, gdzie
ich osobliwe prawa `1/r³` i `1/r⁴` przestają stanowić wiarygodne przybliżenie.

To nadal jest przybliżony model punktowych ładunków: dipole reprezentują
fenomenologicznie spin, a redukcja Landaua–Lifshitza nie zastępuje pełnego
problemu samopola. Stabilność i anihilację prawdziwego pozytonium
opisuje mechanika kwantowa. Program kończy ruch przy `1e-14 m`, ponieważ model
punktowych ładunków nie jest już tam wiarygodny.

## Wymagania i uruchomienie

Wymagany jest CERN ROOT (z `root-config` dostępnym w `PATH`).

```bash
make
```

Przed rozpoczęciem obliczeń program prosi o wybór zjawiska:

```text
1 -> Para-positronium
2 -> Ortho-positronium
3 -> Direct collision
4 -> Scattering
```

W skryptach wybór można przekazać bez interakcji, np. `--phenomenon 2`.

Wykres jest w jednostkach promienia Bohra `a₀`. Odczyt na dole zawiera czas w
pikosekundach, odległość cząstek w pm oraz całkowitą energię mechaniczną układu
elektron–pozyton w eV.
Wyświetlana jest również skumulowana energia wypromieniowana w eV. Linia nad
bieżącymi licznikami pokazuje ich wartości początkowe, a linia pod nimi zmiany
względem początku, oznaczone symbolem `Δ`.

Nad tabelą widoczny jest wynik klasycznego drzewa decyzyjnego oraz cztery
decydujące warunki początkowe: energia względna `E_rel`, orbitalny moment pędu
`L_orb`, przewidywana odległość największego zbliżenia `r_min` i czas życia
`Lifetime`. Dla rozproszenia czas życia jest nieskończony.

Przyciski nad widokiem: `STOP` wstrzymuje animację (zmienia się na `START`), a
`EXIT` natychmiast zamyka program.

Bilans energii cząstek (łącznie z energią dipoli i bliskiego pola Schotta) oraz
zakres odległości można sprawdzić bez uruchamiania okna. Raport nie obejmuje
jeszcze energii pola magnetycznego poruszających się ładunków:

```bash
./positronium --diagnose --phenomenon 1
```

Konkretny losowy przypadek można odtworzyć przez podanie ziarna, np.:

```bash
./positronium --diagnose --phenomenon 4 --seed 8
```

Symulacja jest renderowana w 3D. Przeciągnij myszą w obszarze orbity, aby
zmienić kąt kamery i obejrzeć ruch z innej perspektywy.
