# Semi-klasyczny model pozytonium (CREM-QI)

**CREM-QI** — *Classical Relativistic Electromagnetic Model with Quantum
Imports*: relatywistyczna klasyczna elektrodynamika pary naładowanych cząstek
z jawnie nazwanymi importami kwantowymi (sekcja 5). W kodzie i przełącznikach
model nosi krótką nazwę CREM (`crem_collapse.hpp`, `CREM_*`).

Ten plik opisuje **model fizyczny** w jednym miejscu: z czego się składa, co
jest w nim klasyczne, co jest importem kwantowym i jak wypada wobec pomiarów.
Szczegóły pomiarów, historia zmian i wycofane twierdzenia są w `README.md` i w
pliku audytu `audits/2026-09-10-para-ortho-after-dipole-fix.txt` (z danymi w
`audits/data/section-NNN/`); tutaj są tylko odsyłacze do nich. Stan opisu:
audyt 404 (2026-10-11; 379 dopisany po 380). Każda zmiana modelu trafia tutaj w tym samym
commicie co audyt.

## 1. Czym jest model

Model jest **klasyczną elektrodynamiką dwóch punktowych cząstek**, elektronu i
pozytonu, z których każda niesie ładunek i moment magnetyczny. Do tej
dynamiki dołożono **nazwaną, policzalną listę importów kwantowych**: miejsc,
w których wynik kwantowy wchodzi jako założenie. Model jest semi-klasyczny
dokładnie w tym sensie: trajektoria jest klasyczna, a dyskretność (kwant
działania fotonu, moment pędu fotonu, przygotowanie stanu, reguła anihilacji)
jest wpisana i oznaczona.

Cel jest jeden: sprawdzić, **jak daleko klasyczna elektrodynamika sięga w
odtwarzaniu wyników doświadczalnych pozytonium**, zanim zabraknie jej opisu.
Model ocenia się wobec **pomiarów**, nie wobec innych modeli (sekcja 9).
Żaden człon nie jest dopasowywany do danych; jedynymi wejściami są stałe
CODATA i zmierzone momenty magnetyczne. p-Ps i o-Ps są zawsze liczone i
raportowane osobno.

## 2. Składniki i skale

Każda cząstka ma masę \(m\), ładunek \(\pm e\) i moment magnetyczny
\(\mu=(g/2)\mu_B\) ze zmierzonym \(g\). Para ma masę zredukowaną
\(\mu_{\rm red}=m_e/2\). Wszystkie długości modelu wynikają z jednej skali,
\(\lambda_C\), i jednej liczby bezwymiarowej, \(\alpha\) (relacje są
zablokowane asercjami w `modules/particle_species.hpp`):

| skala | wyrażenie | wartość dla e⁺e⁻ |
|---|---|---|
| sprzężenie Coulomba | \(k=\alpha\hbar c\) | — |
| promień Bohra pary | \(a=a_{\rm pary}=2\lambda_C/\alpha\) | \(2a_0=105{,}8\) pm |
| energia wiązania \(n=1\) | \(B=k/(2a)\) | \(6{,}80285\) eV |
| promień „bariery” | \(r^*=\mu/(ec)=(g/2)\lambda_C/2\) | \(193{,}3\) fm |
| rozmycie pola dipola | \(\varepsilon=0{,}96682\,r^*\) | \(186{,}9\) fm |

\(r^*\) jest dokładnie promieniem, na którym energia dipol–dipol równa się
kulombowskiej: ich stosunek to \((r^*/r)^2\). Długość rozmycia \(\varepsilon\)
nie jest parametrem dopasowania: moment magnetyczny traktuje się jak pętlę
prądu o najmniejszym promieniu, przy którym ładunek \(e\) niesie \(\mu\) z
prędkością nie większą niż \(c\) (\(R=2r^*\)), a profil Plummera ma tę samą
wagę kontaktową co jednorodny dysk tej pętli (audyt 86).

**Zmienne orbity.** Orbita Keplera ma energię \(E=-B/n^2\) i moment pędu
\(L\); \(n=\sqrt{B/|E|}\) jest ciągłe. Główne działanie Delaunaya to
\(J=J_r+L=n\hbar\), a \(e^2=1-(L/n\hbar)^2\). Energia zależy **tylko od
\(J\)**: \(E=-B(\hbar/J)^2\), a częstość orbity to \(\omega=\partial E/\partial J
=2B/(\hbar n^3)\).

## 3. Dynamika klasyczna

Całkowane są położenia, pędy i wektory momentów magnetycznych obu cząstek
(`modules/electrodynamics.hpp`, `modules/crem_engine.hpp`).

1. **Kinematyka relatywistyczna:** \(\mathbf p=\gamma m\mathbf v\); warunek
   \(v<c\) wynika z tej relacji, nie z obcinania.
2. **Oddziaływanie ładunków:** pełne, retardowane pola Liénarda–Wiecherta
   drugiej cząstki i siła Lorentza
   \(\mathbf F_i=q_i[\mathbf E_j+\mathbf v_i\times\mathbf B_j]_{t_{\rm ret}}\),
   dokładne we wszystkich rzędach \(v/c\). Przybliżenie Darwina występuje
   tylko w diagnostyce energii.
3. **Oddziaływanie momentów:** pole dipola rozmyte profilem Plummera
   (\(\varepsilon\) wyżej), z członem magnetyzacji (kontaktowym), tak żeby
   \(\nabla\cdot\mathbf B=0\). Daleko daje to zwykłą energię
   \(U_{dd}\propto[\boldsymbol\mu_1\cdot\boldsymbol\mu_2-3(\boldsymbol\mu_1\cdot\hat r)(\boldsymbol\mu_2\cdot\hat r)]/r^3\),
   której znak dla orbit kołowych zależy od \(P_2(\cos\theta)\). Przy
   \(r\lesssim r^*\) dominuje człon kontaktowy, który **przyciąga momenty
   równoległe w każdym ułożeniu** (audyt 312).
4. **Ruch momentów:** precesja BMT w lokalnym polu. Część orbitalna precesji
   jest identyczna dla obu cząstek, a część od dipola partnera wynosi
   \(\pm1{,}38\cdot10^{11}\) rad/s, ze znakiem zależnym od kanału (audyt 313).
5. **Transport sekularny** (estymator `estimateCremCollapse`, między
   fotonami): wektory \(\mathbf L,\mathbf S_1,\mathbf S_2\) są przenoszone po
   orbicie uśrednionej z zachowaniem \(\mathbf J=\mathbf L+\mathbf S_1+\mathbf S_2\)
   (`modules/secular_spin_orbit.hpp`), z dwiema regułami bilansu:
   - **\(E_{\rm orb}+U\) stałe** (audyt 334): zmiana uśrednionej energii
     dipolowej \(U\) między checkpointami jest odejmowana od energii orbity
     (`CREM_NO_SPIN_ENERGY_EXCHANGE=1` wyłącza). Od audytu 378 \(U\) ma
     jedną gęstość kontaktową z resztą modelu: człon Fermiego
     \(-\tfrac23\mu_0 n_{\rm QR}\,\mathbf m_1\cdot\mathbf m_2\) przy
     \(|L|<\hbar\) i brak kontaktu przy \(l\ge1\); część tensorowa jest
     energią pola uśrednioną po orbicie bez jej części Fermiego–Plummera, a
     przy \(|L|<\hbar\) jest pomijana tak jak w transporcie (izotropia, 357).
     Skok \(U\) wywołany samą emisją zabiera foton: \(E_\gamma=\Delta_{\rm
     drabina}+U_{\rm przed}-U_{\rm po}\) (minus odrzut), a orbita zostaje na
     drabinie (`CREM_NO_PHOTON_SPIN_JUMP=1` oddaje skok orbicie jak w 334,
     `CREM_ORBIT_CONTACT_SPIN=1` przywraca cały stary kontakt);
   - **\(|\mathbf L|\le L_{\rm kol}(E)=n\hbar\)** (audyt 337): nadmiar momentu
     pędu przechodzi do spinów (`CREM_NO_L_CAP=1` wyłącza).
   Sprzężenie \(\mathbf L\cdot\mathbf S\) obraca \(\mathbf L\), ale nie zmienia
   \(|\mathbf L|\); zmienia je tylko część tensorowa spin–spin:
   \(d|L|/dt=-\partial U/\partial g\) (\(g\) – argument perycentrum).
   Przy \(n=1\) część tensorowa daje tylko drgania \(|L|\) o amplitudzie
   \(<0{,}07\hbar\) (audyt 351, rachunek sekularny) — **\(L\) nie miesza się**.
   **Całkowanie (audyty 357–358):** część wspólna precesji wzdłuż \(\mathbf L\)
   (typu \(\mathbf L\cdot\mathbf S\)) jest stosowana jako obrót dokładny
   \(\mathbf S_i(t)=R_{\mathbf J}(A|J|t)\,R_{\mathbf S}(-A|S|t)\,\mathbf S_i(0)\)
   (\(|L|\) stałe), reszta — rozszczepieniem Stranga. Wcześniejszy obrót
   każdego spinu wokół własnej prędkości dawał systematyczny dryf
   \(|L|\propto\theta^3\) (\(\sim0{,}02\)–\(0{,}04\hbar\)/ns przy \(n=1\)); dryfy
   z audytów 351 i 356 były tym błędem (`CREM_SPIN_LEGACY_ROTATION=1`
   przywraca stary obrót, `CREM_SPIN_SUBSTEP` ustawia kąt podkroku).
   **Izotropia stanu \(s\) (audyt 357; domyślna,
   `CREM_NO_S_STATE_ISOTROPY=1` wyłącza):** przy \(|L|<\hbar\) (\(l=0\) w siatce
   Langera) część precesji od dipola partnera,
   \(\boldsymbol\omega_i=\boldsymbol\omega_i^{\rm orb}+\mathsf A_i\mathbf m_j\),
   jest zastąpiona izotropową — stan \(s\) nie ma płaszczyzny, więc tensor
   uśrednia się do zera. Od audytu 378 ta izotropowa część to człon Fermiego
   z gęstością QR, \(\boldsymbol\omega_i=-\gamma_i\tfrac23\mu_0 n_{\rm QR}\mathbf m_j\)
   (\(3{,}7\cdot10^{11}\) rad/s przy \(n=1\)), zamiast śladu
   \(\operatorname{tr}\mathsf A_i/3\) z orbity (kontakt Plummera,
   \(6{,}4\cdot10^{-4}n_{\rm QR}\), audyt 377); przy \(|L|\ge\hbar\)
   uśredniona część Fermiego–Plummera jest odejmowana. Sprzężenie jest typu
   Heisenberga, więc zachowuje \(\mathbf S_1\cdot\mathbf S_2\) i \(w\). Usuwa
   librację kąta wzajemnego w stanie \(s\) (test walidacji `s-state-isotropy`).
6. **Promieniowanie** (sekcja 4).

**Konwencja kanałów.** W p-Ps spiny są antyrównoległe, więc momenty
magnetyczne elektronu i pozytonu są **równoległe**
(\(|\boldsymbol\mu_1+\boldsymbol\mu_2|=2\mu\)). W o-Ps momenty są
antyrównoległe. Kanały różnią się na starcie tylko znakiem
\(\boldsymbol\mu_2\). W statystycznych eksperymentach 1 i 2 spiny są **domyślnie
skwantowane** (\(\cos=\pm1\); audyt 362, decyzja autora; `--free-spins`
przywraca swobodny kąt z audytu 91, przy którym kanał jest tylko
klasyfikacją losowania). Poza nimi domyślny jest kąt swobodny, a
`--spin-quantization` narzuca \(\cos=\pm1\). Waga kanału \(2\gamma\) to
\(w=(|\boldsymbol\mu_1+\boldsymbol\mu_2|/2\mu)^2\). W o-Ps (momenty
antyrównoległe, układ stabilny) \(w=0\) do \(10^{-8}\). W p-Ps klasycznie
(bez transportu singletowego) \(w\) **oscyluje** między \(\approx0{,}62\) a
\(1\) z okresem kilku ps, ze średnią czasową \(\approx0{,}87\)–\(0{,}93\)
(audyty 344, 351, 353); z transportem singletowym (domyślnym przy
skwantowanych spinach, audyt 355) \(w=1\) stale, jak w QM.

**Źródło oscylacji \(w\) (audyty 354–355).** Klasyczna energia spin–spin
\(H=c\,\mathbf S_1\cdot\mathsf T\cdot\mathbf S_2\) z symetrycznym tensorem
\(\mathsf T\) (kontakt + część tensorowa) daje
\(\dot{\mathbf S}_1=c\,\mathbf S_1\times\mathsf T\mathbf S_2\),
\(\dot{\mathbf S}_2=c\,\mathbf S_2\times\mathsf T\mathbf S_1\). Dla
antyrównoległych spinów (\(\mathbf S_1=-\mathbf S_2=\mathbf s\), p-Ps) obie
pochodne są równe, więc \(\mathbf S_1+\mathbf S_2\) rośnie od zera, chyba że
\(\mathsf T\mathbf s\parallel\mathbf s\) (izotropowy kontakt albo
\(\mathbf s\) na osi głównej); dla równoległych (o-Ps) spiny obracają się
sztywnie. Część tensorowa uśredniona po **płaskiej** orbicie jest
anizotropowa: \(|\Delta\omega_\perp|=1{,}2\cdot10^{12}\) rad/s przy
\(L=\hbar/2\) (okres 5,2 ps), \(1{,}5\cdot10^{11}\) rad/s przy \(L=\hbar\);
część orbitalna jest dla obu cząstek identyczna. Przyczyną nie są przeciwne
stosunki giromagnetyczne (\(\gamma_1\gamma_2\) wchodzi tylko do \(c\)), lecz
to, że klasyczna para antyrównoległa ma kierunek \(\mathbf s\), który
\(\mathsf T\) obraca, a kwantowy singlet jest obrotowo niezmienniczy. W QM
każde symetryczne \(S_1\cdot\mathsf T\cdot S_2\) zostawia singlet singletem
(sprawdzone na 200 losowych \(\mathsf T\): odchylenie \(3\cdot10^{-16}\);
niesymetryczne \(\mathsf T\) go psuje), przy każdym \(l\).

**Transport singletowy (audyt 355; domyślny dla p-Ps ze skwantowanymi
spinami, `CREM_NO_SINGLET_TRANSPORT=1` wyłącza):** w transporcie sekularnym
oba momenty singletu obracają się wspólną prędkością
\((\boldsymbol\omega_1+\boldsymbol\omega_2)/2\) — część orbitalna bez zmian,
przeciwne części od partnera się znoszą; \(\mathbf L=\mathbf J-\mathbf S\), więc
\(\mathbf J\) jest zachowane. To import jednego faktu kwantowego (symetryczne
sprzężenie zachowuje \(S=0\)) i **jest konieczny**: klasyczna symetria
zamiany \(P:(\mathbf S_1,\mathbf S_2)\to(\mathbf S_2,\mathbf S_1)\) chroni tylko
swój zbiór punktów stałych \(\mathbf S_1=\mathbf S_2\) (tryplet — dlatego o-Ps
jest stabilne **bez importu**). Zbiór singletowy \(\mathbf S_1=-\mathbf S_2\)
jest punktem stałym \((\mathbf S_1,\mathbf S_2)\to(-\mathbf S_2,-\mathbf S_1)\),
które odwraca nawias Poissona spinów, więc przeprowadza trajektorie w
odwrócone w czasie i niczego nie chroni. W QM \(P\) działa liniowo, a jego
podprzestrzeń antysymetryczna dla dwóch spinów \(\tfrac12\) jest
jednowymiarowa — klasyczne „\(\mathbf S_1=-\mathbf S_2\)” to cała rodzina
kierunków. Pęd ukryty momentu (\(\boldsymbol\mu\times\mathbf E/c^2\)) nie zmienia
tego: wchodzi w sprzężenie spin–orbita, które jest symetryczne (kwantowo
\(\mathbf L\cdot(\mathbf S_1+\mathbf S_2)\), zero na singlecie; człon
\(\mathbf L\cdot(\mathbf S_1-\mathbf S_2)\) znika przy równych masach i \(|g|\)).
Wymiana \(E_{\rm orb}+U\) nie dotyka wektorów spinów. Ograniczenie
\(|\mathbf L|\le L_{\rm kol}\) **zepsułoby** singlet (wkłada nadmiar do
\(\mathbf S_1+\mathbf S_2\)), ale przy transporcie singletowym \(\mathbf J=\mathbf L\)
i nie zadziałało ani razu (audyt 355). Brak klasycznego odpowiednika trypletu
\(m=0\). Efekt: \(w=1{,}000000\) w p-Ps przez
całą kaskadę i trzymanie przy \(n=1\), a \(L\) w p-Ps stoi dokładnie.

**Przygotowanie stanu związanego.** Start na \(a_n=n^2a\), domyślnie
\(n=1\), ostro kołowy: \(L=n\hbar\) z konstrukcji (`--level`).
`CREM_INITIAL_ANGULAR_MOMENTUM=J₀` daje ten sam poziom energii przy
\(L=J_0\,n\hbar\), czyli orbitę o \(e^2=1-J_0^2\) (audyt 301).
`--microcanonical-start` (audyt 326) losuje \(J_0^2\) jednostajnie z \([0,1]\):
to klasyczny zespół mikrokanoniczny przy energii \(E_n\), czyli odpowiednik
całej powłoki \(n\). W zmiennej \(L\) gęstość wynosi
\(p(L)=2L/(n\hbar)^2\) na \([0,n\hbar]\). Strumień losowy jest ten sam co
przy okręgu, więc oba starty da się porównywać parami.

**Start Langera** (audyt 353; przez `CREM_INITIAL_ANGULAR_MOMENTUM` z
\(J_0=(l+\tfrac12)/n\)): \(L=(l+\tfrac12)\hbar\), \(J_r=(n_r+\tfrac12)\hbar\),
\(n=n_r+l+1\) — kwantowanie półklasyczne WKB z indeksem Maslowa (import).
Reguły fotonu zmieniają \(L\) o \(\pm\hbar\) albo odbijają je
(\(L\to\hbar-L\)), więc część ułamkowa \(f=(L/\hbar)\bmod1\) przechodzi
\(f\to f\) lub \(f\to1-f\): sieć \(L\in\{\tfrac12,\tfrac32,\dots\}\hbar\) jest
**zachowywana**, a \(f=\tfrac12\) jest punktem stałym odbicia. Każda para
startująca z siatki Langera dochodzi do \(n=1\) z \(L=\hbar/2\)
(zmierzone: \(|L-\tfrac12|\le8\cdot10^{-3}\)). Dynamika tego \(L\) **nie
wybiera** (brak atraktora, audyt 351) — wybiera je przygotowanie.

## 4. Promieniowanie: dwa tryby

**Ciągły:** zredukowana siła Landaua–Lifshitza na każdą cząstkę
(`--radiation-reaction individual`); domyślny w trybie wizualnym i w
wiązkach. Orbita traci energię płynnie po obwiedni Larmora; nie ma fotonów
ani kwantu, a orbita spada aż do bariery (klasyczna katastrofa promienista).

**Fotonowy (`--radiation-reaction stochastic`; od audytu 341 domyślny w
statystycznych eksperymentach 1 i 2, czyli dla stanów związanych):** moc
E1 trafia do hazardu i jest wypłacana dyskretnymi fotonami, które zachowują
pęd (estymator sekularny `estimateCremCollapse`, `modules/crem_collapse.hpp`).
Reguły, w kolejności działania:

- **Hazard** (wyłącznie fotonowy od audytu 328; strata mierzonej orbity też
  trafia do hazardu, `CREM_CONTINUOUS_ORBIT_CREDIT=1` przywraca kredyt
  ciągły): \(\lambda=S(e)\,P_{E1}(a,e)/E_{\rm ref}\), gdzie \(E_{\rm ref}\) to
  energia fotonu z reguły niżej, a \(S(e)\) uwzględnia rozkład mocy orbity
  eliptycznej na harmoniczne (tablica interpolowana w \(e^2\), przy małym
  \(e\): \(S\approx1-2e^2\); audyt 340). Domyślnie foton pada, gdy skumulowany
  hazard osiągnie 1 (deterministycznie); `--emission poisson` losuje próg z
  \(\mathrm{Exp}(1)\). Moc \(P_{E1}\) liczy się po węzłach orbity z siłami
  opóźnionymi; równanie Keplera rozwiązuje Newton zabezpieczony przedziałem
  \(E-M\in[-e,e]\) (audyt 348; goły Newton rozbiegał się przy \(e\ge0{,}985\)).
- **Energia — reguła działania (domyślna od audytu 352;
  `CREM_NO_ACTION_PHOTON=1` przywraca \(k\hbar\omega\)).** Foton zabiera
  \(\hbar\) głównego działania \(J\). Ponieważ \(E\) zależy tylko od \(J\),
  \[\Delta E_{\rm orb}=\int_{J-\hbar}^{J}\omega\,dJ=E(n)-E(n-1)\]
  **przy dowolnym mimośrodzie**. Harmoniczna \(k\), losowana z widma dipola
  orbity Keplera (`eccentricOrbitHarmonicNumber`), wybiera
  \(\Delta n=\min(k,\lfloor n-1\rfloor)\) — zawsze całkowita liczba kwantów
  działania (audyt 355; wcześniej podłoga \(n'=1\) dopuszczała ułamkowe
  \(\Delta J\), gdy wymiana energii spinowej odsunęła \(n\) od całkowitego):
  \(\Delta=B[1/(n-\Delta n)^2-1/n^2]\), tolerancja \(10^{-4}\) na \(n\). Dla
  całkowitych \(n\) identyczne z poprzednią regułą.
  Odrzut: orbita traci dokładnie \(\Delta\), więc foton dostaje
  \(E_\gamma=\Delta-\Delta^2/(2W)\), gdzie \(W=Mc^2+\mu_{\rm red}E\).
  Klasyczne \(E_\gamma=\hbar\omega\) jest granicą tej reguły dla małego skoku.
- **Kierunek i skrętność:** kierunek z wzoru kątowego pola dalekiego,
  \(\propto1+\cos^2\theta\) względem \(\mathbf L\); skrętność \(h=\pm1\) z
  \(P(+)=(1+c)^2/(2(1+c^2))\).
- **Moment pędu (\(\Delta l=-1\), od audytu 330):** foton zabiera \(\hbar\)
  wzdłuż osi orbity, \(\mathbf L'=\mathbf L-\hbar\hat{\mathbf L}\), czyli
  \(|L'|=\big||L|-\hbar\big|\) (`CREM_DIRECTIONAL_PHOTON_SPIN=1` przywraca
  \(\hbar h\) wzdłuż kierunku lotu). Przy \(|L|\ge\hbar\) daje to
  \(\Delta J_r=0\) (działanie radialne zachowane, jak w kwantowej kaskadzie
  yrast). Przy \(|L|<\hbar\) reguła **odbija** \(L\) do \(\hbar-|L|\) — tego
  QM nie ma (z \(l=0\) dozwolone jest tylko \(\Delta l=+1\)); wariant testowy
  `CREM_DL_PLUS_BELOW_HBAR=1` daje wtedy \(L'=L+\hbar\) (audyt 348–349).
- **Sufit energii:** orbita z \(L'\) musi istnieć (\(e^2\ge0\)), czyli
  \(E'\ge-B(\hbar/L')^2\). Przy regule działania foton, którego sufit nie
  mieści, jest **odrzucany**, a nie przycinany (przycięcie złamałoby
  \(\Delta J=\hbar\)); przy \(n\le1\) każdy foton jest odrzucany.
- **Między fotonami:** dynamika zachowawcza i transport sekularny (sekcja 3).
- **Stan końcowy (`EmissionChannelClosed`):** przy regule działania kanał
  jest zamknięty, gdy nie ma całego kwantu (\(\lfloor n-1\rfloor=0\), czyli
  \(n<2\)) albo gdy najmniejszy krok (\(k=1\)) nie zmieści \(L'\):
  \(L'>(n-1)\,\hbar\). Bez reguły działania warunek
  brzmi \(|L/\hbar-1|\ge n\) (audyt 308).

**Co daje reguła działania** (audyt 349, start mikrokanoniczny \(n=2\),
oba kanały identycznie):
- każdy foton \(2\to1\) ma \(E_\gamma=5{,}10212\) eV, a każda kaskada ląduje
  dokładnie na \(n=1\) z \(L'=|L_0-\hbar|\);
- **poziomy \(n<1\) nie istnieją**: z \(n=1\) nie da się zabrać \(\hbar\)
  działania, więc nie ma rodzin kołowych \(n\approx0{,}44\) ani eliptycznych
  z audytów 343–345, które były artefaktami reguły \(k\hbar\omega\);
- **start przy `--level 1` nie emituje nic** (stan końcowy od razu);
  kaskadę daje dopiero \(n\ge2\). Kaskada \(2\to1\) w programie głównym
  (start mikrokanoniczny, 4 pary) trwa średnio \(9{,}2\pm2{,}7\) ns w obu
  kanałach i potrzebuje `--crem-wallclock-budget-s` rzędu 1200 s — przy
  domyślnych 90 s wszystkie przebiegi są cenzurowane (audyt 352);
- dyskretność jest **siecią niezmienniczą, nie atraktorem**: całkowite \(n\)
  zostaje całkowite, ale model nie wybiera całkowitych \(n\) sam — robi to
  przygotowanie stanu;
- z \(\Delta l=+1\) przy \(|L|<\hbar\) stany \(L_0<\hbar\) przy \(n=2\) są
  zamknięte (odpowiednik metastabilnego \(2s\), bo E1 \(2s\to1s\) jest
  zabronione) i \(J_r\) jest zachowane do \(5\cdot10^{-5}\).

## 5. Importy kwantowe

Każdy import ma przełącznik albo jest nazwany jako nieusuwalny. Ich
historyczne ceny są w README („Sześć importów kwantowych i ich cena”).

Pełna lista (stan audytu 369), z oceną, czy import da się uzyskać z
dynamiki modelu. Wejścia mierzone (\(\alpha\), masy, \(\mu\), \(g\)) nie są
importami.

| # | import | co wpisuje | domyślnie | czy z dynamiki? |
|---|---|---|---|---|
| 1 | \(\hbar\) działania na foton (reguła działania) | \(\Delta J=-\hbar\): skoki \(E(n)-E(n-\Delta n)\), \(n=1\) końcowe | tak (`CREM_NO_ACTION_PHOTON`) | nie — to sam kwant |
| 2 | \(\hbar\) momentu pędu fotonu | \(\mathbf L'=\mathbf L-\hbar\hat{\mathbf L}\) (\(\Delta l=-1\)) | tak | **nie wynika z #1** poza \(e=0\): klasycznie na \(\hbar\) działania przypada \(f(e)\hbar\) momentu pędu; schemat odwrotny (tempo z działania, \(\Delta l=\pm1\) z momentu siły) odrzucają czasy życia wodoru (audyt 370) |
| 3 | \(\Delta l=+1\) przy \(|L|<\hbar\) | reguła wyboru z \(l=0\) | tak (`CREM_DL_REFLECT_BELOW_HBAR`) | nie |
| 4 | tempo fotonów \(f(e)P/\hbar\omega=\langle\dot L\rangle_{\rm klas}/\hbar\) | liczba fotonów na czas | tak (audyty 368, 370; `CREM_HAZARD_SPECTRAL_TABLE` — tablica \(S(e)\)) | **tak — wynika z #2** (moment siły promieniowania w jednostkach \(\hbar\)); wodór 2p / 3p: +5,1% / +1,7% od pomiaru |
| 5 | całkowite \(n\) na starcie | \(a_n=n^2a\) (`--level`) | tak | nie (reguła działania zachowuje całkowite \(n\), nie wybiera go) |
| 6 | start Langera | \(L=(l+\tfrac12)\hbar\) | tak w eksp. 1/2 (`--circular-start`) | zachowywany przez dynamikę (punkt stały odbicia), ale nie wybierany; jedyne przypisanie z dokładnym \(\langle r^{-2}\rangle\) (kontakt 1s i \(\langle p^4\rangle\), niżej); cena: \(\langle r^{-3}\rangle\) \(-11{,}1\%\) przy \(l=1\) niezależnie od \(n\) (spin–orbita 2P, audyt 386); „\(+\tfrac12\)” to faza Maslova — efekt falowy |
| 7 | kwantyzacja orientacji spinów | singlet / tryplet \(m=\pm1\) | tak w eksp. 1/2 (`--free-spins`) | nie; długość \(|S_i|=\hbar/2\) wynika ze zmierzonych \(\mu\), \(g\) |
| 8 | transport singletowy | wspólna precesja w singlecie | tak (`CREM_NO_SINGLET_TRANSPORT`) | **nie — udowodnione** (audyt 355g: klasyczna symetria chroni tylko tryplet) |
| 9 | izotropia stanu \(s\) | tensor uśredniony przy \(|L|<\hbar\) | tak (`CREM_NO_S_STATE_ISOTROPY`) | **nie**: płaszczyzna stanu \(s\) się nie obraca (p-Ps: \(\mathbf J=\mathbf L\); o-Ps \(1^3S_1\): \(\mathbf S\parallel\mathbf L\)); **ograniczona** (audyt 389): wartość = ścisła średnia QM dla \(l=0\), próg nie jest parametrem, a bez reguły poziom \(1^3S_1\) zależałby od orientacji o 131,7 GHz (\(\tfrac94\) różnicy kontaktowej) — ~700× ponad zgodność pomiarów HFS w polu i w zerowym polu |
| 10 | tożsamość Quigga–Rosnera | \(n_{\rm kontakt}=\hbar/(2\pi a^3n^3L)\) | tak (`CREM_ORBIT_CONTACT_DENSITY`) | nie (twierdzenie QM; klasyczna jest tylko \(\langle dV/dr\rangle\)) |
| 11 | \(\sigma v=4\pi r_e^2c\) | tempo anihilacji przy kontakcie | tak | nie (QED) |
| 12 | \(\varepsilon_{\rm OP}\) | tłumienie \(3\gamma\) | tak | nie (QED) |
| 13 | rozkład płaszczyzny \(3\gamma\) | \(1-\tfrac13\cos^2\theta_n\) | tylko eksperyment 6 | nie (QED) |
| 14 | ZPF \(\hbar\omega/2\) | pole punktu zerowego (SED) | nie (opcja) | — ; nie stabilizuje \(n=1\): bez stanu stacjonarnego, jonizacja przez pompowanie mimośrodu (audyt 387) |

**Typy fotonów (audyt 396).** Każdy foton ma typ i to typ decyduje, co
niesie. E1 (j = 1, Δl = ±1) pochodzi z dipola orbity, M1 (j = 1, Δl = 0)
z precesji momentów własnych, E2 (j = 2, Δl = 0, ±2, s → s zabronione)
z kwadrupola ładunkowego orbity, κ₂ = (q₁m₂² + q₂m₁²)/M², zerowego dla e⁺e⁻.
Kanał losowany jest z udziału tempa; E2 ma własne widmo harmonicznych i rzut
m z FFT kwadrupola orbity Keplera (suma = postać zamknięta typu Petersa
P_E2 = κ₂²k³f(e)/(5πε₀c⁵μ³a⁵)), rozkład kątowy j = 2 i przekaz mħ orbitalnego
L. M1 (audyt 397) nie rusza orbity: odwraca jeden moment (|ΔS| = ħ) i niesie
U_przed − U_po sprzężenia spinów; w stanach s to dokładnie HFS modelu
(4/3)μ₀μ²n_QR, a gdy żadne odwrócenie nie obniża U (p-Ps w stanie s), foton
jest odmawiany. E3 i M2 (β⁴) są pominięte. Dla pozytonium
liczby się nie zmieniają (κ₂ = 0, udział M1 ~10⁻¹⁵); w wodorze udział E2 jest
rzędu α² (6,5·10⁻⁶ przy 3d) i daje m.in. bezpośrednie 3d → 1s.

**Jedno założenie skalowe.** Skala całego modelu bierze się z \(\hbar\):
przy ustalonym \(L=\hbar\) minimum \(L^2/(2\mu r^2)-k/r\) leży dokładnie w
\(a\) (argument Bohra), a reguła działania przenosi to samo \(\hbar\) na
skoki energii. Stan 1s ma kwantowo \(l=0\), a jego skalę daje energia
lokalizacji, której model nie ma.

**Czy całkowite \(n\) da się wyprowadzić? Zegar wewnętrzny (audyt 380).**
Moment magnetyczny w modelu to pętla prądowa o promieniu \(R=2r^*=\hbar/(mc)\),
po której ładunek obiega z prędkością \(c\). To zegar o częstości
\(\Omega=mc^2/\hbar\), czyli zegar de Broglie'a. Na orbitach Keplera
pozytonium daje on trzy wyniki:
- zegary \(e^-\) i \(e^+\) biegną dokładnie synchronicznie (\(v_1=v_2\) w
  układzie środka masy), więc nie ma między nimi dudnień;
- poślizg dylatacji czasu na orbitę wynosi w rzędzie wiodącym \(J/(4\hbar)\)
  cykli na zegar, przy każdej ekscentryczności (wirial: \(\langle v^2\rangle\)
  nie zależy od \(e\)). Zegar liczy więc **główne działanie \(J\)**, a nie
  \(L\) — tę samą zmienną, którą skwantowuje reguła działania. Zależność od
  \(e\) pojawia się dopiero w następnym rzędzie,
  \(\delta=(\alpha^2/16n^2)\,(4J/L-3)\) (z \(\langle v^4\rangle\)). Wyjaśnia ona
  zmierzone odchyłki do \(3\cdot10^{-5}\) z resztą \(\le5\cdot10^{-9}\), a po
  dodaniu wyrazu \(v^6\) — \(\le2\cdot10^{-9}\) (audyt 380e). Czynnik
  \(4J/L-3\) ma postać relatywistycznej poprawki kinetycznej
  Sommerfelda \(n/l-\tfrac34\), bo faza zegara to działanie swobodnej
  cząstki \(-mc^2\!\int d\tau/\hbar\);
- faza opóźnienia \(\Omega r/c\) między zegarami nie ma żadnej struktury przy
  całkowitym \(L/\hbar\).

Zablokowanie fazy zegara na tle o stałej częstości dałoby \(J=4N\hbar\), a
nie \(J=n\hbar\). Warunek Bohra wymaga fazy \(2\pi J/\hbar\) na orbitę, czyli
fali o liczbie falowej \(\mu v/\hbar\) — fali fazowej de Broglie'a, której
model nie ma. Reguła działania pozostaje importem. Brakującym elementem jest
dokładnie ta fala: zegar dostarcza 1/4 fazy (1/2 dla obu zegarów), a fala
dostarczyłaby resztę.

**Spójność startu Langera (import 6) z importem 10 i poprawką \(\langle p^4\rangle\).**
Przy \(L=(l+\tfrac12)\hbar\) i \(J=n\hbar\) klasyczna średnia po orbicie Keplera
\(\langle r^{-2}\rangle=1/(n^3(l+\tfrac12)a^2)\) jest **dokładnie** równa
kwantowej. Z tej jednej tożsamości wynikają dwa trafienia:
- gęstość kontaktowa QR (\(\propto\langle dV/dr\rangle\propto\langle r^{-2}\rangle\))
  daje ściśle \(|\psi_{ns}(0)|^2\) (audyty 350, 361);
- relatywistyczna poprawka kinetyczna \(\langle p^4\rangle/(8m^3c^2)\) ma
  klasycznie czynnik \(4J/L-3\), czyli ściśle kwantowe \(4n/(l+\tfrac12)-3\)
  w \(\Delta E=-(E_n^2/2mc^2)[4n/(l+\tfrac12)-3]\). Ten sam czynnik daje drugi
  rząd dylatacji zegara (audyt 380e).

Żadnego importu to nie usuwa, bo oba trafienia są skutkiem importów 5 i 6, a
nie ich zamiennikiem. Import 10 też zostaje: związek „gęstość w zerze =
\((\mu/2\pi\hbar^2)\langle dV/dr\rangle\)” jest twierdzeniem kwantowym,
ponieważ klasyczna orbita nie ma gęstości w zerze.

Daje to natomiast uzasadnienie wartości „\(+\tfrac12\)”. Langer to jedyne
przypisanie \(L\), przy którym \(\langle r^{-2}\rangle\) jest dokładne, więc
jednocześnie trafia w kontakt \(1s\) i w \(\langle p^4\rangle\).
\(\sqrt{l(l+1)}\) i \(l\) nie trafiają w żadne z nich. Nie jest to jednak
wyprowadzenie, a pomiary nie wybierają jednego przypisania: czasy życia
wodoru 2p–5p różnią się dla Langera o \(+5\) do \(+20\%\), a
\(\sqrt{l(l+1)}\) wypada lepiej przy \(n=4,5\) (audyty 371–372). Samo
„\(+\tfrac12\)” to faza Maslova (\(\pi/2\) na punkt zwrotny), czyli efekt
falowy wymagający fali o liczbie falowej \(p/\hbar\). Tej fali model nie ma
(audyt 380).

**Cena Langera: tabela średnich (audyt 386).** Każda średnia po orbicie
wybiera inne \(L\), więc jedno przypisanie trafia tylko w część z nich:

| średnia | błąd przy \(L=(l+\tfrac12)\hbar\) | \(L\), które ją uściśla | gdzie w modelu |
|---|---|---|---|
| \(\langle r^{-1}\rangle\) | 0 (wirial) | każde | energia |
| \(\langle r^{-2}\rangle\) | 0 | \((l+\tfrac12)\hbar\) | kontakt QR (anihilacja, Fermi HFS, ZPF), \(\langle p^4\rangle\) |
| \(\langle r^{-3}\rangle\) | \(l(l+1)/(l+\tfrac12)^2-1\): \(-11{,}1\) / \(-4{,}0\) / \(-2{,}0\%\) przy \(l=1,2,3\), bez zależności od \(n\) | \(\sqrt[3]{l(l+\tfrac12)(l+1)}\,\hbar\) | spin–orbita, tensor spin–spin |
| \(\langle r\rangle\) | 1s \(-8{,}3\%\), 2p \(-2{,}5\%\) | \(\sqrt{l(l+1)}\,\hbar\) | — |
| \(\langle r^2\rangle\) | 1s \(-29\%\), 2s \(-8{,}3\%\) | \(\sqrt{l(l+1)-\tfrac13}\,\hbar\) (brak dla \(l=0\)) | — |

Stała spin–orbita 2P w silniku (2,437 GHz) to dokładnie
\(\tfrac89\cdot\)Breit (2,7416 GHz) z dokładnością \(3\cdot10^{-5}\). Z
pomiarów \(2^3P\) wychodzi 2,743 GHz, więc cały deficyt \(-11{,}1\%\) to
czynnik Langera. Zgodność odstępu \(P_2-P_1\) do 1% (audyt 374) jest
kompensacją dwóch błędów: spin–orbity \(-0{,}61\) GHz i za słabego
klasycznego tensora \(+0{,}56\) GHz. Dla stanów \(ns\) żaden rozkład \(L\)
przy stałym \(a\) nie daje kwantowego \(\langle r^2\rangle\): klasyczne
maksimum wynosi \(5n^4/2\) wobec \((n^2/2)(5n^2+1)\), czyli 0,83 przy
\(n=1\). Dla 1s kwantowe \(\langle r\rangle\) wymaga \(L=0\), a wtedy
kontakt jest nieskończony. Rozmiar stanu s wymaga rozrzutu energii. Pytanie
„czy dynamika wybiera \(L=\hbar/2\)” jest więc przy stałym \(a\) źle
postawione. Test dynamiczny musi dotyczyć rozkładu \((E,L)\).

**Dopasowanie faz: skąd całkowite \(\Delta n\), a skąd nie całkowite \(n\).**
Wszystkie klasyczne mechanizmy promieniowania ładunku spełniają jeden warunek
dopasowania faz, \(\omega-\mathbf k\cdot\mathbf v=s\,\omega_{\rm orb}\):
- promieniowanie Czerenkowa to \(s=0\) w ośrodku, przy \(v>c/n_r\);
- undulator i cyklotron to \(s\ge1\).

Foton pary **nie** jest promieniowaniem Czerenkowa, z trzech powodów. Para
jest w próżni (\(n_r=1\)) i porusza się z \(v\approx\alpha c\). Promieniuje
ładunek przyspieszany, a nie ładunek w ruchu jednostajnym: tempo zgadza się z
Larmorem do \(10^{-4}\) (audyt 371). Rozkład kątowy jest dipolowy
\((1+\cos^2\theta)\), bez stożka.

Dla orbity związanej warunek sprowadza się do \(\omega=s\,\omega_{\rm orb}\):
foton jest \(s\)-tą harmoniczną ruchu okresowego. **Liczba \(s\) jest
całkowita dzięki okresowości ruchu, a nie kwantyzacji.** Reguła działania
wiąże ją z liczbą przeskoczonych szczebli: harmoniczna \(k\) zabiera
\(k\hbar\), czyli \(\Delta n=k\) (zasada korespondencji,
`actionStepQuanta`). Klasyczna fizyka daje więc całkowite **\(\Delta n\)**,
ale nie całkowite **\(n\)**. To, od którego szczebla zaczyna się drabina, czyli
dlaczego \(J=n\hbar\), wymaga fazy \(2\pi J/\hbar\) na orbitę. Taką fazę daje
fala de Broglie'a, a nie dopasowanie faz promieniowania (audyt 380).

**Całkowite \(l\).** W modelu całkowite \(l\) powstaje tak samo jak
całkowite \(n\): jest wpisane na starcie i zachowywane, ale nie wybierane.
- Start Langera (import 6) ustawia \(L=(l+\tfrac12)\hbar\) z całkowitym \(l\):
  \(l=n-1\) albo `--orbital-l`.
- Każdy foton zmienia \(L\) o dokładnie \(\hbar\) wzdłuż osi orbity: import 2
  (\(\Delta l=-1\)) i, poniżej \(\hbar\), import 3 (\(\Delta l=+1\)). Siatka
  \((l+\tfrac12)\hbar\) jest więc zachowana. Start z \(L=0{,}8\hbar\) dałby
  siatkę \(0{,}8\to1{,}8\to\dots\)
- **\(l\le n-1\) wynika klasycznie**, bez importu: na orbicie Keplera
  \(L\le J\) (koło ma największe \(L\) przy danej energii), więc przy
  \(L=(l+\tfrac12)\hbar\), \(J=n\hbar\) mamy \(l+\tfrac12\le n\).
- \(l\ge0\) to minimum siatki, \(L=\hbar/2\), którego pilnuje import 3.

W QM całkowite \(l\) bierze się z jednoznaczności fali na sferze, czyli z
warunków EBK na kątach: \(\oint p_\varphi d\varphi=2\pi m\hbar\) (bez punktów
zwrotnych) i \(\oint p_\theta d\theta=2\pi(l-|m|+\tfrac12)\hbar\) (dwa punkty
zwrotne, faza Maslova). Razem dają \(L=(l+\tfrac12)\hbar\), czyli przepis
Langera. To ten sam warunek fali stojącej co dla \(n\), tylko na kącie, więc
wymaga tej samej fali \(p/\hbar\), której model nie ma.

**Wspólny wzór dla \(n\) i \(l\): klasycznie całkowite są różnice, a nie
wartości.**
- \(\Delta n\): foton to \(s\)-ta harmoniczna ruchu okresowego w czasie
  (\(s\) całkowite), a reguła działania wiąże \(s\) z \(\Delta n\).
- \(\Delta l\): klasyczne pole multipolowe o liczbie azymutalnej \(m\) niesie
  moment pędu i energię w stosunku \(m/\omega\) (Jackson, rozdz. 9.8). \(m\)
  jest całkowite z jednoznaczności pola w \(\varphi\), więc przy porcji
  energii \(\hbar\omega\) porcja momentu pędu wynosi \(m\hbar\); dla E1
  \(\Delta l=\pm1\). Ściśle tylko na kole: przy \(e>0\) porcja działania
  \(\hbar\) nie oznacza \(\hbar\) momentu pędu (audyt 370), dlatego import 2
  jest osobnym importem.

Wartości \(n\) i \(l\), czyli punkt startu siatki, wymagają fazy fali
stojącej: \(2\pi J/\hbar\) na orbitę i \((l+\tfrac12)2\pi\) na kąt. Model je
wpisuje na starcie (importy 5 i 6).

## 6. Anihilacja

**Klasyczna dynamika nie zawiera anihilacji.** Model ma dwie jawnie
kwantowe recepty.

### 6.1. Tempo w estymatorze (audyt 342; tryb fotonowy)

\[\Gamma=\sigma v\,[\,w+(1-w)\,\varepsilon_{\rm OP}\,]\,\langle n\rangle_{\rm orb},\qquad
E[T]=\int_0^\infty S(t)\,dt,\quad S=e^{-\int\Gamma dt}\]

- \(\sigma v=4\pi r_e^2c\) (import; z \(|\psi_{1s}(0)|^2=1/(\pi a^3)\) daje
  \(T_0=1/(\sigma v|\psi_{1s}(0)|^2)=124{,}49\) ps);
- \(\varepsilon_{\rm OP}=4(\pi^2-9)\alpha/(9\pi)=1/1113{,}9\) — tłumienie
  \(3\gamma\) (Ore–Powell, wiodący rząd QED);
- \(\langle n\rangle_{\rm orb}\) — gęstość kontaktowa uśredniona po orbicie z
  jądrem Plummera \(3\varepsilon^2/(4\pi\rho^5)\).

**Ta recepta nie pasuje do stanu \(n=1\)** (audyty 343–347, 349d): orbita
\(n=1\) z \(L\sim\hbar\) prawie nie dochodzi do zera, więc
\(\langle n\rangle_{\rm orb}\) wynosi \(10^{-6}\)–\(10^{-3}|\psi_{1s}(0)|^2\), a
czasy życia wychodzą \(10^4\)–\(10^5\) razy za długie. \(|\psi_{1s}(0)|^2\)
odtwarza tylko wąskie \(L^*=0{,}178\hbar\), które jednak przesuwa się z
\(\varepsilon\) jak \(\varepsilon^{0{,}28}\), czyli jest dopasowaniem, a nie
przewidywaniem (audyt 347). W języku reguły działania kontakt to \(J=0\)
(„\(n=0\)”), a anihilacja to ostatni krok \(1\to0\), którego tempa estymator
jeszcze nie ma.

### 6.2. Krok \(1\to0\) z twierdzenia Quigga–Rosnera (audyt 350; w silniku od audytu 361)

**W silniku (audyt 361; `CREM_ORBIT_CONTACT_DENSITY=1` przywraca 6.1):**
gęstość kontaktowa w tempie anihilacji to \(n_{\rm QR}\) przy \(L<\hbar\) i
\(0\) przy \(l\ge1\) (\(L\ge\hbar\)); program główny drukuje w podsumowaniu
eksperymentów 1/2 czas życia w stanie końcowym \(1/\Gamma\), rozrzut
\(\langle1/\Gamma\rangle\langle\Gamma\rangle\) i \(E[T]\) od przygotowania.
Wymaga skwantowanych spinów — **domyślnych w eksperymentach 1/2 od audytu
362**; ze swobodnym kątem (`--free-spins`) \(w\ne0,1\) i o-Ps wychodzi
\(\sim0{,}5\) ns (audyt 361). Program główny bez flag daje p-Ps 124,494 ps i
o-Ps 138,674 ns, oba jednowykładniczo (start \(1s\)); ze startem \(2p\)
(`--level 2`) te same czasy w stanie końcowym, a \(E[T]\) od przygotowania
większe o kaskadę \(\approx10\) ns.

**Kinematyka:** fotony anihilacyjne (\(2\gamma\)/\(3\gamma\)) zabierają cały
\(\mathbf J=\mathbf L+\mathbf S_1+\mathbf S_2\), więc reguła \(\Delta l=-1\)
tego kroku nie ogranicza; kanał wybiera \(w\).

**Tempo:** twierdzenie Quigga–Rosnera (QM, ścisłe dla stanów s):
\(|\psi(0)|^2=(\mu_{\rm red}/2\pi\hbar^2)\langle dV/dr\rangle\). Klasyczna
orbita ma \(\langle k/r^2\rangle=k\hbar/(n^3a^2L)\), więc
\[n_{\rm QR}=\frac{1}{2\pi a^3n^3}\,\frac{\hbar}{L},\qquad
\tau(L)=\frac{1}{\Gamma}=T_0\,\frac{2n^3L/\hbar}{w+(1-w)\varepsilon_{\rm OP}}.\]
Przy \(n=1\) i \(L=\hbar/2\) (Langer) \(n_{\rm QR}=|\psi_{1s}(0)|^2\). Wynik
nie zależy od \(\varepsilon\); importem jest tożsamość QR, nie
\(|\psi(0)|^2\). Przy \(n=2\) czynnik \(n^3=8\) daje
\(\tau_{2s}=8\tau_{1s}\) dla tego samego \(L\), tak jak
\(|\psi_{2s}(0)|^2=|\psi_{1s}(0)|^2/8\).

**Średnie po zespole.** Start mikrokanoniczny przy \(n=2\): w jednostkach
\(\hbar\) \(L_0=2j_0\), \(j_0^2\sim U(0,1)\), czyli \(p(L_0)=L_0/2\) na
\([0,2]\). Średni czas życia mieszaniny to
\(\langle\tau\rangle=\int p(L)\,\tau(L)\,dL\propto\langle L'\rangle\).

- **A′** (odbicie, \(L'=|L_0-1|\), wszystkie pary na \(n=1\)):
  \[\langle L'\rangle=\int_0^2|x-1|\,\frac{x}{2}\,dx
  =\int_0^1\frac{(1-x)x}{2}dx+\int_1^2\frac{(x-1)x}{2}dx=\frac1{12}+\frac5{12}=\frac12,\]
  więc \(\langle\tau\rangle_p=T_0=124{,}5\) ps (przy \(w=1\)) i
  \(\langle\tau\rangle_o=T_0\cdot1113{,}9=138{,}7\) ns.
- **C′** (\(\Delta l=+1\) przy \(L_0<1\)): zespół dzieli się na dwie
  populacje. \(P(L_0<1)=\int_0^1x/2\,dx=1/4\) zostaje przy \(n=2\) (\(2s\));
  \(3/4\) schodzi na \(n=1\) z \(L'=L_0-1\):
  \[\langle L'\rangle=\frac{\int_1^2(x-1)x\,dx}{\int_1^2x\,dx}=\frac{5/6}{3/2}=\frac59
  \;\Rightarrow\;138{,}3\ \text{ps},\ 154\ \text{ns}.\]
  Populacja \(2s\) ma \(\langle L_0\rangle=\int_0^1x^2dx/\int_0^1x\,dx=2/3\),
  więc \(\langle\tau\rangle_{2s}=T_0\cdot16\cdot\tfrac23=10{,}7\,T_0\).

**Wagi.** A′ i C′ to **alternatywne reguły**, nie populacje — nie uśrednia
się ich. Wewnątrz wariantu wagi zespołu startowego są już w
\(\langle L'\rangle\); w C′ jeden czas życia całego zespołu wymagałby wag
\(3/4\) (\(n=1\)) i \(1/4\) (\(n=2\)). Pomiary \(125{,}14\) ps i \(142{,}04\)
ns dotyczą stanu podstawowego, więc porównuje się populację \(n=1\). p-Ps i
o-Ps się nie waży: to osobne składowe widma czasów życia.

**Zbieżność, nie wyprowadzenie.** \(\langle L'\rangle=1/2\) w A′ jest
wartością Langera, ale wynika z dwóch arbitralnych wyborów (zespół
mikrokanoniczny przy \(n=2\), reguła odbicia). Mieszanina **nie jest
jednowykładnicza**: gęstość \(L'\) przy zerze jest skończona, więc
\(\langle1/L'\rangle\) (tempo początkowe) rozbiega się logarytmicznie, a ogon
\(S(t)\) ma czas 1,5–1,9× średniej. Zmierzone rozpady są jednowykładnicze,
a \(L\) nie miesza się w czasie życia (audyt 351).

**Zmierzone na trajektoriach** (start mikrokanoniczny \(n=2\), 48 par, stany
\(n=1\) bez bariery):

| wariant | silnik | p-Ps \(\langle\tau\rangle\) | o-Ps \(\langle\tau\rangle\) | \(\tau_o/\tau_p\) | \(\langle1/\Gamma\rangle\langle\Gamma\rangle\) p / o |
|---|---|---|---|---|---|
| A′ | 349/350 (bez reguł spinowych) | 160,8 ps (1,28×; \(w\approx0{,}92\)) | 150,0 ns (1,06×) | 933 | 1,60 / 1,67 |
| C′ | 349/350 | 187,5 ps (1,50×) | 177,6 ns (1,25×) | 947 | 1,41 / 1,36 |
| A′ | **obecny (audyt 359)** | **137,4 ps** (1,10×; \(w=1\), \(\langle L\rangle=0{,}552\)) | **150,1 ns** (1,06×) | 1092 (z barierą 1114) | 1,52 / 1,67 |
| C′ | **obecny (audyt 359)** | **159,4 ps** (1,27×) | **177,6 ns** (1,25×) | 1114 | 1,36 / 1,36 |

Przy \(w=1\) średnia p-Ps to dokładnie \(T_0\cdot2\langle L\rangle\)
(\(124{,}49\times2\times0{,}552=137{,}4\) ps), a stosunek kanałów —
\(1/\varepsilon_{\rm OP}\), o ile obie próby mają ten sam skład (w A′ bez
bariery próby różnią się o 1–2 trajektorie, stąd 1092). Rozpad pozostaje
niejednowykładniczy, bo rozrzut \(L\) pochodzi ze startu.

**Ze startem Langera** (audyt 353; \(2p\), \(3p\), \(3d\) i \(2s\) w A′;
26 par na kanał na \(n=1\), wszystkie z \(L=\hbar/2\)):

| kanał | \(\langle\tau\rangle\) | mediana | \(\langle1/\Gamma\rangle\langle\Gamma\rangle\) | pomiar |
|---|---|---|---|---|
| p-Ps | 136,3 ps | 127,65 ps | 1,0225 (oscylacja \(w\)) | 125,14 ps |
| o-Ps | 138,65 ns | 138,67 ns | 1,00002 (jednowykładniczo) | 142,04 ns |

**Z transportem singletowym** (audyt 355, \(2p\), 6 par): p-Ps \(w=0{,}999996\)
w każdej parze, \(\tau=124{,}495\) ps w każdej (\(-0{,}52\%\) od pomiaru,
\(\langle1/\Gamma\rangle\langle\Gamma\rangle=1{,}000000\)); o-Ps bez zmian (138,66 ns,
wiersze identyczne); \(\tau_o/\tau_p=1113{,}78\) (pomiar 1135, \(-1{,}9\%\)).

Bez tej reguły: o-Ps różni się od pomiaru o \(-2{,}4\%\), czyli tyle, ile
\(1/\varepsilon_{\rm OP}=1113{,}9\) od \(1135\). Rozrzut p-Ps pochodzi
wyłącznie z oscylacji \(w\) (min 0,59, średnia 0,934), nie z \(L\). Stan
\(2s\) w C′ (zamknięty przy \(n=2\)) ma \(\tau_{\rm QR}=995{,}95\) ps (p-Ps)
i \(1{,}109\) µs (o-Ps), czyli \(8\times\) więcej, jak
\(|\psi_{2s}(0)|^2=|\psi_{1s}(0)|^2/8\).

### 6.3. Eksperyment 6 (zdarzenie anihilacji; od audytu 375 w stanie końcowym)

`--mode statistical --phenomenon 6 --runs N`: para jest przygotowana jak w
eksperymentach 1/2 (Langer \(L=(l+\tfrac12)\hbar\), \(l=n-1\), domyślnie
\(n=1\), więc \(L=\hbar/2\); spiny skwantowane, transport singletowy dla
p-Ps). Silnik prowadzi ją do stanu końcowego (kaskada fotonów, reguła
działania) i podaje to samo tempo, które całkują eksperymenty 1/2:
\[
\Gamma=\sigma v\,[w+(1-w)\varepsilon_{\rm OP}]\,n_{\rm QR},\qquad
t=t_{\rm kaskady}+\mathrm{Exp}(1/\Gamma).
\]
Czas \(t\) jest losowany z własnego strumienia pary, od audytu 379 z przeżycia \(S(t)\), które silnik całkuje wzdłuż kaskady (odcinkami stałe \(\Gamma\)), a za stanem końcowym z \(\mathrm{Exp}(1/\Gamma)\) (`drawDecayFromHazard`; sekcja 6.4). Stan s powyżej \(n=1\) może więc anihilować w trakcie kaskady. Kanał wynika z
\(P(2\gamma)=w/(w+(1-w)\varepsilon_{\rm OP})\), a fotony powstają w układzie
spoczynkowym pary. Dla \(2\gamma\) oś jest izotropowa. Dla \(3\gamma\) energie
mają dokładny rozkład Ore–Powella, a płaszczyzna rozpadu spełnia
\(dN/d\cos\theta_n\propto1-\tfrac13\cos^2\theta_n\) względem
\(\mathbf S_1+\mathbf S_2\) (drzewowa amplituda QED). Oś spinu stanu
końcowego jest izotropowa: \(|S|=0\) dla p-Ps, \(\hbar\) dla o-Ps.
Na wyjściu są: rozkład czasów rozpadu w każdym kanale, \(\tau\) (MLE),
test KS wobec rozkładu wykładniczego, stosunek \(\tau_{\rm o}/\tau_{\rm p}\)
i fotony. Eksperyment nie wnosi nowych importów; korzysta z importów 6.1–6.2
i 13 z sekcji 5.

Audyt 375 (1000 par) daje:

| | p-Ps | o-Ps |
|---|---|---|
| \(1/\Gamma\) z silnika | 124,494 ps | 138,674 ns |
| \(\tau\) z próby | \(118{,}9\pm3{,}8\) ps | \(148{,}3\pm4{,}7\) ns |
| test KS | \(p=0{,}90\) | \(p=0{,}99\) |
| kanał | 100% \(2\gamma\) | 100% \(3\gamma\) |

Stosunek z próby wynosi \(1248\pm4{,}5\%\) i jest zgodny z \(1/\varepsilon=1113{,}9\)
w granicach szumu (pomiar: 1135,0). Fizyczna treść to wartości silnika; próba
jest ich realizacją zdarzeniową.

Stary mechanizm „przy barierze” (audyty 316–320) działa pod
`--contact-barrier` i daje wyjście identyczne bajt w bajt. Anihilacja następuje
tam przy pierwszym wejściu w \(r\le r^*\) przy dynamice zachowawczej i starcie
\(J_0=0{,}07\). Czasy bezwzględne nie mają w nim sensu, a stosunek \(\approx2\cdot10^4\)
zawiera klasyczny czynnik geometrii kontaktu bez odpowiednika kwantowego.

**Energia fotonów anihilacji (audyt 379).** Fotony dzielą niezmiennik
\(W=(m_1+m_2)c^2-B_{\rm k}+U\), gdzie \(U\) to ta sama energia sprzężenia
spinów co w bilansie (378): \(\mp1{,}21\cdot10^{-4}\) eV dla p-Ps / o-Ps
przy \(n=1\). Wcześniej była tu energia dipolowa w perycentrum,
\(\sim\pm10^{-2}\) eV, ze znakiem zależnym od orientacji (człon tensorowy,
który w stanie s usuwa reguła izotropii).

### 6.4. Definicja czasu życia pary

**Czas życia pary** to czas \(T\) od przygotowania pary (chwila \(t=0\),
stan startowy) do jej anihilacji, czyli zniknięcia pary i powstania fotonów
\(2\gamma\) lub \(3\gamma\). \(T\) jest zmienną losową. Jej rozkład wyznacza
tempo anihilacji \(\Gamma(t)\) z sekcji 6.2, liczone na trajektorii, którą
silnik prowadzi przez kaskadę fotonów:
\[
\Gamma(t)=\sigma v\,\bigl[w(t)+(1-w(t))\,\varepsilon_{\rm OP}\bigr]\,n_{\rm QR}(t),\qquad
S(t)=P(T>t)=\exp\Bigl(-\int_0^t\Gamma\,dt'\Bigr),\qquad
f(t)=\Gamma(t)\,S(t).
\]
- \(n_{\rm QR}=\hbar/(2\pi a_B^3n^3L)\) dla \(|L|<\hbar\) i 0 dla \(l\ge1\).
  Stany p (i wyższe \(l\)) nie anihilują; anihilują stany s, przy każdym
  \(n\).
- \(w(t)\) to waga \(2\gamma\) z momentów magnetycznych w danej chwili. Przy
  spinach skwantowanych jest stała: \(w=1\) dla p-Ps, \(w=0\) dla o-Ps.
- W silniku \(\Gamma\) jest **odcinkami stałe**: jedna wartość na checkpoint,
  przeliczana po każdym fotonie i zmianie orbity
  (`annihilationHazardSegments`). Za ostatnim checkpointem obowiązuje tempo
  stanu końcowego \(\Gamma_{\rm k}\).
- Kanał rozpadu losuje się w chwili \(T\) z
  \(P(2\gamma)=w/(w+(1-w)\varepsilon_{\rm OP})\).

**Zegar.** Wszystkie czasy są czasem własnym pary, mierzonym w jej
chwilowym układzie spoczynkowym (`simulatedTimeTotal`). Odpowiedniki
laboratoryjne (`lifetimeSecondsLab`) różnią się o czynnik \(\gamma\) ruchu
środka masy po odrzutach fotonów. Przy prędkości odrzutu \(\sim10^{-5}c\) to
różnica \(\sim10^{-11}\) względnie.

**Wielkości wypisywane i ich znaczenie** (każdy kanał osobno; trajektorie
p-Ps i o-Ps są sparowane po ziarnie):

| wielkość | definicja | gdzie | z czym porównywać |
|---|---|---|---|
| **czas życia stanu końcowego** \(\tau_{\rm k}\) | \(1/\Gamma_{\rm k}\), tempo w stanie, na którym kaskada się kończy (\(n=1\), \(L=\hbar/2\)); w wydruku średnia \(\langle1/\Gamma_{\rm k}\rangle\) po trajektoriach, z rozrzutem \(\langle1/\Gamma\rangle\langle\Gamma\rangle\) (1 = jeden wykładnik w całym zespole) | eksperymenty 1/2, „lifetime at the final state” | **zmierzone czasy życia** (125,14 ps, Al-Ramadhan i Gidley 1994; 142,04 ns, Vallery i in. 2003) |
| **średni czas od przygotowania** \(E[T]\) | \(\int_0^\infty S\,dt\): całka po odcinkach kaskady plus \(S(t_{\rm k})/\Gamma_{\rm k}\) | eksperymenty 1/2, „E[T] from the preparation” | wielkość modelu; wliczony czas kaskady |
| **czas kaskady** \(t_{\rm k}\) | czas od przygotowania do zamknięcia kanału emisji (stan końcowy) | `lifetimeSeconds`, histogramy „collapse time” | to **nie** jest czas życia pary: para po \(t_{\rm k}\) nadal istnieje |
| **czas rozpadu jednej pary** \(t\) | losowanie z \(S(t)\): z odcinków kaskady, a za stanem końcowym z \(\mathrm{Exp}(1/\Gamma_{\rm k})\) (audyt 379; w 375 \(t=t_{\rm k}+\mathrm{Exp}(1/\Gamma_{\rm k})\)) | eksperyment 6 | \(\tau\) (MLE) = średnia \(t\) w kanale = estymator \(E[T]\); rozkład \(t\) i test KS |

**Związki między nimi.**
- Przy starcie Langera z \(n=1\) kaskady nie ma: \(t_{\rm k}=0\), \(S(t_{\rm k})=1\)
  i \(E[T]=\tau_{\rm k}\). Rozkład \(T\) jest wtedy jednym wykładnikiem.
- Przy starcie ze stanu p (\(l=n-1\)) kaskada przechodzi przez stany bez
  kontaktu, więc \(E[T]\approx t_{\rm k}+\tau_{\rm k}\). Na przykład p-Ps z
  \(n=2\): \(3{,}49\) ns \(=3{,}36+0{,}12\) ns.
- Przy starcie ze stanu s powyżej \(n=1\) para może anihilować przed końcem
  kaskady. Zmierzone dla 3s (audyt 381): p-Ps anihiluje w stanie 3s w 39 z 40 par,
  średnio po \(2{,}6\pm0{,}4\) ns, czyli \(\approx1/\Gamma_{3s}=3{,}4\) ns. Tempo z \(2s\) jest \(8\times\), a z \(3s\) \(27\times\) mniejsze niż
  z \(1s\). \(E[T]\) i rozkład \(T\) to uwzględniają, a \(\tau_{\rm k}\) nie.

**Dlaczego z pomiarem porównuje się \(\tau_{\rm k}\), a nie \(E[T]\).**
Zmierzony czas życia to stała zaniku wykładniczego ogona widma czasowego, a
nie średni czas od powstania pary. Ten ogon pochodzi od pozytonium w stanie
podstawowym, więc jego odpowiednikiem w modelu jest \(1/\Gamma\) stanu
\(1s\).

**Czym czas życia nie jest.**
- Nie jest **czasem zapadania** klasycznej orbity do bariery Comptona ani do
  granicy retardacji. Tak nazywano wynik trybu z ciągłą reakcją promienistą
  (`--radiation-reaction individual`) przed audytem 342. Tamta wielkość
  nadal nazywa się w kodzie „collapse time” i opisuje utratę energii orbity,
  a nie anihilację.
- Nie jest **czasem wejścia w \(r\le r^*\)** ze starego eksperymentu 6
  (`--contact-barrier`, audyty 316–320). Tam czasy bezwzględne nie mają sensu
  fizycznego.
- Dotyczy pary w próżni. Anihilacja na elektronach ośrodka (pick-off) i
  gaszenie polem nie są modelowane; gaszenie polem 50 µT jest pomijalne
  (audyt 376).

**Trajektorie niedokończone.** Trajektoria, która nie dojdzie do stanu
końcowego w budżecie czasu rzeczywistego (`--crem-wallclock-budget-s`;
domyślnie 90 s, a 1200 s przy starcie z \(n\ge2\)), jest **cenzurowana**. Nie
daje \(\tau_{\rm k}\) ani \(E[T]\), a wydruk podaje ich liczbę. Histogramy
czasu kaskady traktują ją jako obserwację cenzurowaną prawostronnie
(Kaplan–Meier) i niczego nie ekstrapolują.

## 7. Ruch par p-Ps i o-Ps — podsumowanie badań (audyty 320–359)

Ta sekcja opisuje, co para robi w modelu w obecnych ustawieniach domyślnych
(tryb fotonowy, reguła działania, spiny skwantowane), krok po kroku i
**osobno dla kanałów**, oraz jak wypada to wobec pomiarów. Szczegóły i
wyprowadzenia są w sekcjach 3–6.

### 7.1. Stan początkowy

- **Orbita:** elipsa Keplera na poziomie \(n\) (całkowite; \(n=1\) jest
  stanem końcowym, więc kaskada wymaga \(n\ge2\)), półoś \(a_n=n^2a\).
  Moment pędu: **start Langera** \(L=(l+\tfrac12)\hbar\) (jedno \(L\) na
  \(l\); import WKB, audyt 353; **domyślny w eksperymentach 1/2 od audytu
  361**, \(l=n-1\), `--orbital-l`, `--circular-start`) albo zespół mikrokanoniczny
  \(p(L)=2L/(n\hbar)^2\) (rozrzut \(L\)). Obie cząstki krążą wokół wspólnego
  środka masy po przeciwnych stronach, z prędkością względną
  \(v\approx\alpha c/n\).
- **Spiny:** p-Ps — \(\mathbf S_1=-\mathbf S_2\) (momenty magnetyczne
  **równoległe**, \(w=1\)); o-Ps — \(\mathbf S_1=\mathbf S_2\) (momenty
  **antyrównoległe**, \(w=0\)). o-Ps odpowiada \(m=\pm1\); trypletu \(m=0\)
  klasycznie nie ma. **Orbita startowa jest w obu kanałach taka sama**;
  kanały różnią się tylko znakiem \(\boldsymbol\mu_2\).

### 7.2. Ruch orbitalny między fotonami

- Ruch względny jest zachowawczy (pełne siły Lorentza, opóźnione pola,
  siły dipolowe); w estymatorze sekularnym elementy orbity \(n\), \(|L|\) są
  między fotonami stałe — poza wymianą energii z oddziaływaniem dipolowym
  (\(E_{\rm orb}+U\) stałe, przesuwa \(n\) o \(\sim10^{-5}\)).
- **Płaszczyzna orbity i perycentrum precesują:** sprzężenie
  spin–orbita obraca \(\mathbf L\) i \(\mathbf S\) wokół \(\mathbf J\) przy
  stałym \(|L|\) (obrót dokładny, audyt 358); perycentrum precesuje w
  płaszczyźnie. \(\mathbf J=\mathbf L+\mathbf S_1+\mathbf S_2\) jest zachowane
  do błędu zaokrągleń.
- **\(|L|\) się nie zmienia.** Część tensorowa spin–spin daje przy \(l\ge1\)
  tylko drgania \(<0{,}07\hbar\); przy \(l=0\) (\(L<\hbar\)) usuwa ją reguła
  izotropii stanu \(s\) (audyt 357). Wcześniej obserwowane dryfy \(L\) (351,
  356) były błędem całkowania (audyty 357–358).

### 7.3. Ruch spinów

| | p-Ps (singlet) | o-Ps (tryplet \(m=\pm1\)) |
|---|---|---|
| precesja orbitalna (spin–orbita) | wspólna dla obu spinów | wspólna dla obu spinów |
| pole dipola partnera, klasycznie | rozsuwa spiny (anizotropia części tensorowej): \(w\) oscyluje 0,6–1, okres 5 ps przy \(L=\hbar/2\) (354) | obraca oba spiny sztywnie, \(w=0\) dokładnie — **bez importu** (symetria zamiany chroni \(\mathbf S_1=\mathbf S_2\)) |
| w modelu (domyślnie) | **transport singletowy** (import, 355): wspólna średnia prędkość, \(w=1{,}000000\) | bez zmian, \(w=0\) |
| stan \(s\) (\(L<\hbar\)) | tylko izotropowa część pola partnera (357) | j.w. |

### 7.4. Emisja fotonów — kaskada

- **Tempo (audyty 368, 370; `CREM_HAZARD_FROM_GAP=1` przywraca 352,
  `CREM_HAZARD_SPECTRAL_TABLE=1` tablicę \(S(e)\)):** liczba fotonów na
  jednostkę czasu to klasyczne tempo utraty momentu pędu w jednostkach
  \(\hbar\), \(f(e)\,P_{E1}/\hbar\omega_{\rm orb}\) z \(f(e)=\omega\,dL/dE
  =(1-e^2)^{3/2}/(1+e^2/2)\), a każdy
  foton niesie przerwę z reguły działania. Klasyczna moc Larmora jest zegarem
  emisji, nie źródłem energii (energia uchodzi w tempie
  \(P\cdot E_{\rm przerwa}/\hbar\omega\)). Daje to zmierzone czasy życia wodoru:
  2p 1,681 ns (pomiar 1,600(4)), 3p 5,675 ns (5,58(13)) przy liniach 10,199 i
  12,087 eV (audyty 367–370).
- **Każdy foton zabiera \(\hbar\) działania** (reguła działania, 352):
  \(\Delta E_{\rm orb}=E(n)-E(n-\Delta n)\), \(\Delta n=\min(k,\lfloor
  n-1\rfloor)\), oraz \(\hbar\) momentu pędu wzdłuż osi orbity
  (\(\Delta l=-1\)); przy \(L\ge\hbar\) działanie radialne \(J_r\) zostaje
  zachowane, jak w kaskadzie yrast.
- **Linie** (z odrzutem): \(2\to1\) \(5{,}10212\) eV w każdej kaskadzie,
  \(3\to1\) \(6{,}04698\) eV, \(3\to2\) \(0{,}94484\) eV; \(3d\to2p\to1s\)
  (dwa fotony), \(3p\to1s\) bezpośrednio (353).
- **\(|L|<\hbar\)** (\(l=0\)): wariant C′ (**domyślny od audytu 361**,
  \(\Delta l=+1\)) zostawia \(2s\) **metastabilne**, zgodnie z QM; wariant A′
  (`CREM_DL_REFLECT_BELOW_HBAR=1`) odbija \(L\to\hbar-L\), przez co
  \(2s\to1s\) zachodzi przez E1 — niezgodnie z QM (348–349, 353).
- **Czas kaskady:** \(2p\to1s\) **3,36 ns** (audyt 370; 3,33 ns w 368) (wcześniej
  \(\approx10\) ns z hazardem dzielonym przez przerwę; teoria najniższego rzędu
  dla Ps: 3,187 ns, pomiaru brak — Day i in. 1992 bez wyniku); dla wodoru
  model daje 2p 1,681 ns i 3p 5,675 ns wobec zmierzonych 1,600 i 5,58 ns.
  **Taki sam w obu kanałach** — fotony E1 nie zależą od spinów
  (w C′ ta sama przyczyna stopu i liczba fotonów w 48/48 par).

### 7.5. Stan końcowy \(n=1\)

- Brak całego kwantu działania: żaden foton nie wychodzi, \(n=1\) jest stanem
  końcowym (stany \(n<1\) nie istnieją).
- \(L=\hbar/2\) przy starcie Langera (sieć \(L=(l+\tfrac12)\hbar\) jest
  zachowywana przez reguły fotonu, \(f=\tfrac12\) to punkt stały odbicia;
  zmierzone \(|L-\tfrac12|\le1\cdot10^{-4}\)); przy starcie mikrokanonicznym
  \(L=|L_0-\hbar|\) jest rozrzucone.
- \(L\) i \(w\) pozostają stałe w czasie życia (351, 355, 357–358).

### 7.6. Anihilacja — krok \(1\to0\)

Tempo z tożsamości Quigga–Rosnera (sekcja 6.2):
\(\tau=T_0\cdot2n^3(L/\hbar)/[w+(1-w)\varepsilon_{\rm OP}]\). Przy
\(L=\hbar/2\): p-Ps (\(w=1\)) \(\tau=T_0\); o-Ps (\(w=0\))
\(\tau=T_0/\varepsilon_{\rm OP}\). Kanał wybiera \(w\): p-Ps \(\to2\gamma\),
o-Ps \(\to3\gamma\) (tłumienie \(\varepsilon_{\rm OP}\)).

### 7.7. Zgodność z pomiarami

| obserwabla | p-Ps: model | p-Ps: pomiar | o-Ps: model | o-Ps: pomiar | ocena |
|---|---|---|---|---|---|
| czas życia, start Langera (2p, audyt 358; program główny bez flag, audyt 362) | 124,494 ps | 125,14 ps | 138,674 ns | 142,04 ns | \(-0{,}52\%\) / \(-2{,}37\%\) |
| kształt rozpadu, start Langera | jednowykładniczy | jednowykładniczy | jednowykładniczy | jednowykładniczy | zgodny |
| czas życia, start mikrokanoniczny (A′, 359) | 137,4 ps | 125,14 ps | 150,1 ns | 142,04 ns | \(+9{,}8\%\) / \(+5{,}6\%\); kształt **niezgodny** (niejednowykładniczy) |
| \(\tau_o/\tau_p\) (wspólny dla kanałów) | 1113,97 | 1135,0 | — | — | \(-1{,}9\%\) = \(1/\varepsilon_{\rm OP}\) (QED wiodącego rzędu, import) |
| linia \(2\to1\) (\(2p\to1s\), audyt 378) | 5,102244 eV | — | 5,102001 eV | 5,10179 eV (1S–2S, tryplet) | o-Ps \(+0{,}21\) meV (\(2P\), nie \(2S\): porównanie przybliżone o strukturę subtelną \(n=2\)); różnica kanałów 0,243 meV = HFS 1s modelu (58,53 GHz); przed 378 obie 5,10212 eV |
| kanał rozpadu | \(2\gamma\) (\(w=1\)) | \(2\gamma\) | \(3\gamma\) (\(w=0\)) | \(3\gamma\) | zgodny (reguła \(\varepsilon_{\rm OP}\) importowana) |
| struktura subtelna \(n=2\), \(\nu_J=E(2^3S_1)-E(2^3P_J)\) | \(\nu_0,\nu_1,\nu_2=-6{,}88,\ -9{,}86,\ -14{,}19\) GHz | \(+18{,}50,\ +13{,}01,\ +8{,}62\) GHz (Hagena i in. 1993) | — | — | model kładzie \(2^3S_1\) **poniżej** \(2^3P_J\) (różnica 22,8–25,4 GHz); brak dwóch członów stanów S bez klasycznego odpowiednika (wirtualna anihilacja 10,95 GHz, kontaktowy Darwin 5,47 GHz) i części spinowej P (~6–9 GHz); odstęp \(P_2-P_1\) zgodny (4,33 wobec 4,39 GHz), \(P_1-P_0\) o 46% za mały (audyt 363) |
| HFS 1s (o-Ps − p-Ps) | \(58{,}53\) GHz | \(203{,}39\) GHz | (różnica kanałów) | | \(0{,}288\) pomiaru, znak poprawny (audyt 360): klasyczne korelacje spinów dają połowę kwantowej części magnetycznej (117,06 GHz), a 3/7 HFS (wirtualna anihilacja) nie ma klasycznego odpowiednika |

Pomiary: \(\tau_p\) — Al-Ramadhan i Gidley, PRL 72, 1632 (1994);
\(\tau_o\) — Vallery i in., PRL 90, 203402 (2003); 1S–2S — Fee i in., PRL 70,
1397 (1993); HFS — Ishida i in., PLB 734, 338 (2014).

### 7.8. Co jest wyprowadzone, a co importowane

**Z modelu (klasyczna elektrodynamika + ustalone reguły):** ruch orbitalny i
spinów, odrzut, zachowanie \(\mathbf J\) i \(J_r\), skala linii przy danym
\(\hbar\), stabilność trypletu o-Ps, zachowanie sieci Langera przez kaskadę,
czynnik \(n^3\) dla \(2s\), równość kaskad w obu kanałach.

**Importy (każdy nazwany, z wyłącznikiem):** \(\hbar\) działania na foton
(352), \(\hbar\) momentu pędu fotonu (330), start Langera (353), kwantyzacja
spinów, transport singletowy (355), izotropia stanu \(s\) (357), tożsamość
Quigga–Rosnera (350), \(\sigma v=4\pi r_e^2c\), \(\varepsilon_{\rm OP}\).

**Wniosek.** Zgodność czasów życia do 0,5–2,4% i kształtu rozpadu opiera się
na tych importach; własnym wkładem modelu jest to, że klasyczna dynamika je
**utrzymuje** (siatka \(L\), \(w\), \(J\), działanie) przez całą kaskadę i
czas życia, a nie że je wyprowadza. Różnice \(-0{,}5\%\) / \(-2{,}4\%\) i
\(+0{,}33\) meV są rzędu poprawek QED wyższego rzędu, których model nie ma.

## 8. Przełączniki testowe i pomiarowe (domyślnie wyłączone)

| przełącznik | działanie | audyt |
|---|---|---|
| `CREM_NO_ACTION_PHOTON=1` | energia fotonu \(k\hbar\omega\) zamiast reguły działania | 352 |
| `CREM_DL_REFLECT_BELOW_HBAR=1` | odbicie \(L\to\hbar-L\) przy \(|L|<\hbar\) (A′) zamiast \(\Delta l=+1\) | 361 |
| `CREM_ORBIT_CONTACT_DENSITY=1` | gęstość kontaktowa z orbity (Plummer) zamiast Quigga–Rosnera | 361 |
| `--circular-start`, `--orbital-l <l>` | start \(L=n\hbar\) zamiast Langera; wybór \(l\) w starcie Langera | 361 |
| `CREM_NO_SINGLET_TRANSPORT=1` | p-Ps bez transportu singletowego (niezależna precesja) | 355 |
| `CREM_NO_S_STATE_ISOTROPY=1` | pełne (anizotropowe) pole partnera także przy \(|L|<\hbar\) | 357 |
| `CREM_SPIN_SUBSTEP=<rad>` | maks. kąt podkroku transportu spin–orbita (domyślnie 0,05) | 357 |
| `CREM_SPIN_LEGACY_ROTATION=1` | stary obrót spinów (bez dokładnej części \(\mathbf L\cdot\mathbf S\); dryf \(|L|\propto\theta^3\)) | 358 |
| `CREM_CLOSE_BELOW_HBAR=1` | kanał zamknięty przy \(|L|<\hbar\) | 346 |
| `CREM_DL_PLUS_BELOW_HBAR=1` | (od 361 bez działania — \(\Delta l=+1\) jest domyślne) | 348–349 |
| `CREM_HOLD_AFTER_CLOSURE=<s>` | po zamknięciu trzymaj parę bez emisji, drukuj \(L\), \(n\), \(w\) | 351 |
| `CREM_NO_L_CAP=1`, `CREM_NO_SPIN_ENERGY_EXCHANGE=1` | wyłączają reguły bilansu z sekcji 3 | 337, 334 |
| `CREM_DIRECTIONAL_PHOTON_SPIN=1` | \(\hbar h\) wzdłuż kierunku fotonu zamiast osi | 330 |
| `CREM_CONTINUOUS_ORBIT_CREDIT=1` | kredyt ciągły straty mierzonej orbity | 328 |
| `--contact-barrier` | eksperyment 6: stary mechanizm anihilacji przy wejściu w \(r^*\) (audyty 316–320) zamiast tempa QR w stanie końcowym | 375 |
| `CREM_ORBIT_CONTACT_SPIN=1` | kontakt Fermiego w transporcie spinów i w \(U\) z orbity (Plummer) zamiast QR; wyłącza też skok energii spinów w fotonie | 378 |
| `CREM_REFUSE_ACTION_STEP=1` | foton obniżający \(L\) przy zbyt dużym \(\Delta n\) odrzucany zamiast skróconego skoku | 384 |
| `CREM_HARMONIC_TABLE=1` | harmoniczna z tablicy kwantyli zamiast dokładnych udziałów \(w_k\) (Bessel) | 384 |
| `CREM_SPIN_S_STATE_EXACT=1` | stany s: przepływ spinów \(R_J R_S\) (L·S + Heisenberg dokładnie, reszta połówkami kroku) zamiast podkroków 0,05 rad — **niedokładny** (2,9e-4 wobec 1e-6 w jednym kroku), tylko do badań | 385 |
| `CREM_SPIN_SLOW_RATES=1` | stany s z QR: dwie średnie BMT na ocenę tempa zamiast jednej (wynik ten sam, wolniej) | 388 |
| `CREM_S_STATE_THRESHOLD=<x>` | próg reguły izotropii w \(\hbar\) (domyślnie 1); test: 0,6 i 1,4 dają wynik bit w bit | 389 |
| `--ps-source rest\|thermal[:T]\|beam:<E>` | ruch środka masy Ps w laboratorium dla fotonów anihilacji i zegara laboratoryjnego eksp. 1/2 (dynamika bez zmian) | 390 |
| `CREM_CONTINUOUS_MUTUAL_RADIATION=1` | tryb `stochastic`: ciągła wymiana radiacyjna przez wzajemne pola opóźnione (sprzed audytu 395; część interferencyjna liczona podwójnie) | 395 |
| `CREM_FORCE_E2=<p>` | wymusza kanał E2 z prawdopodobieństwem p (pomiar; naturalny udział ~α²) | 396 |
| `CREM_DEBUG_MULTIPOLE=1` | wydruk typu, k, m i L każdego fotonu ścieżki sekularnej | 396 |
| `CREM_M1_ORBITAL_QUANTUM=1` | foton M1 z kwantem orbitalnym ħω_orb zamiast odwrócenia spinu (sprzed audytu 397) | 397 |
| `CREM_MECHANICAL_HAZARD_SCALE=<x>` | mnoży hazard fotonów ścieżki mechanicznej (pomiar reguły na orbicie związanej) | 398 |
| `CREM_DEBUG_LADDER=1` | wydruk poziomu i energii każdego fotonu ścieżki mechanicznej | 398 (399: bilans okna) |
| `CREM_WINDOW_SCALE_L=1` | spłata okna przez skalowanie pędu (L odchodzi z energią, dE/dL = ω) zamiast ΔJ, ΔL fotonu | 399 |
| `CREM_LADDER_WINDOW_FRACTION=<ε>` | długość okna drabiny: zmiana J o ε n na obieg (domyślnie 0,01) | 399 |
| `CREM_DEBUG_LADDER_SPLIT=1` | zmiana energii w otwartym oknie drabiny rozbita na przyrost fotonu, obrót L i dynamikę (z członami) | 400 |
| `CREM_LADDER_CLAMP_ROTATION=1` | obrót L w oknie przycina cel do r\|p⊥\| (zeruje p_r) zamiast zostawić pęd | 401 |
| `CREM_MECHANICAL_FORCE_E2=1` | każdy foton ładunkowy ścieżki mechanicznej pary z kwadrupolem jest E2 (pomiar) | 401 |
| `CREM_MECHANICAL_STEPS_PER_ORBIT=<N>` | gęstość kroku ścieżki mechanicznej (domyślnie 128) | 401 |
| `CREM_MECHANICAL_TOLERANCE=<tol>` | tolerancja względna integratora ścieżki mechanicznej (domyślnie 1e-5) | 401 |
| `CREM_VELOCITY_CORRECTOR=1` | siły końca kroku liczone ponownie przy prowizorycznej prędkości końcowej (usuwa błąd rzędu 1 członów zależnych od prędkości; test) | 402 |
| `CREM_INTERACTION_SAMPLE_ONLY=1` | eksp. 5: tylko losowanie (indeks, K, b) bez całkowania — do odnalezienia zdarzenia | 403 |
| `CREM_INTERACTION_EVENT_INDEX=<i>` | eksp. 5: całkuje tylko zdarzenie i (pozostałe zostają NumericalFailure; diagnoza) | 403 |
| `CREM_M1_TORQUE_REDUCED=1` | moment siły reakcji M1 z m‴ = Ω×(Ω×(Ω×m)) (precesja BMT) zamiast z różnic po historii (test) | 404 |
| `CREM_DUMP_FAIL=<plik>` | zrzut stanu i historii pierwszego kroku, który nie spełnił tolerancji (diagnoza offline) | 404 |
| `CREM_LAB_DIPOLE_DRIFT=1` | zegar laboratoryjny z \(\gamma\) zsumowanego dryfu dipolowego (audyt 110) zamiast samego odrzutu fotonów | 382 |
| `CREM_STOP_BELOW_SURVIVAL=<S>` | kończy trajektorię, gdy przeżycie anihilacyjne spadnie poniżej \(S\) (przyczyna zatrzymania `SurvivalThreshold`); do długożyjących stanów s | 381 |
| `CREM_EXP6_PARA_ONLY=1` | eksperyment 6 bez trajektorii o-Ps | 381 |
| `CREM_PERIAPSIS_DIPOLE_ENERGY=1` | energia dzielona przez fotony anihilacji z energią dipolową w perycentrum (stara recepta) zamiast \(U\) z bilansu | 379 |
| `CREM_NO_PHOTON_SPIN_JUMP=1` | skok \(U\) przy emisji idzie do orbity (334), nie do fotonu | 378 |
| `CREM_EMISSION_REACH=1` | wydruk każdego fotonu (\(E\), \(L\), \(e^2\), \(k\), \(n\)) | — |

Pełna lista z uzasadnieniami: `modules/configuration_panel.hpp`.

## 9. Wobec pomiarów

Pomiary odniesienia: \(\tau_{\rm para}=125{,}14\) ps (Al-Ramadhan i Gidley,
PRL 72, 1632 (1994)); \(\tau_{\rm orto}=142{,}04\) ns (Vallery i in., PRL 90,
203402 (2003)) — średnie czasy życia, nie półokresy; odstęp 1S–2S
\(5{,}10179\) eV (Fee i in., PRL 70, 1397 (1993)); HFS \(203{,}3942\) GHz
(Ishida i in., PLB 734, 338 (2014)); \(m_ec^2=510{,}99895\) keV (CODATA 2018).

| obserwabla | model | pomiar | ocena |
|---|---|---|---|
| linia \(2\to1\) | o-Ps \(5{,}102001\) eV, p-Ps \(5{,}102244\) eV w 100% kaskad (reguła działania + skok energii spinów, audyt 378; przedtem \(5{,}10212\) eV w obu) | \(5{,}10179\) eV (1S–2S, tryplet) | o-Ps \(+0{,}21\) meV (\(4\cdot10^{-5}\)); foton \(2p\to1s\) wobec linii \(1S\)–\(2S\), więc porównanie przybliżone o strukturę subtelną \(n=2\); reszta — poprawki subtelne i QED, których drabina \(-B/n^2\) nie ma (audyt 349) |
| stan końcowy kaskady | dokładnie \(n=1\) | \(1s\) | tak; przy starcie mikrokanonicznym \(L=|L_0-\hbar|\) rozrzucone, przy starcie Langera \(L=\hbar/2\) w każdej parze |
| linie \(3\to1\), \(3\to2\) | \(6{,}04698\), \(0{,}94484\) eV (kaskada yrast \(3d\to2p\to1s\), audyt 353) | — | pomiaru nie ma w repozytorium |
| \(\tau\) z kroku \(1\to0\) (QR, start Langera, transport singletowy) | p-Ps 124,495 ps, o-Ps 138,68 ns (audyt 358), oba jednowykładniczo | 125,14 ps / 142,04 ns | \(-0{,}52\%\) / \(-2{,}4\%\) (audyt 355); bez transportu singletowego p-Ps średnia \(+8{,}9\%\) |
| \(2s\) metastabilne | w C′ brak fotonu z \(2s\) (E1 zabronione) | QM: tak | zgodne; w A′ \(2s\) schodzi przez odbicie (niezgodne z QM) |
| \(\tau\) z recepty 6.1 przy \(n=1\) | \(10^4\)–\(10^5\times\) za długie | 125,14 ps / 142,04 ns | recepta nie pasuje do \(n=1\) (sekcja 6.1) |
| \(\tau\) z kroku \(1\to0\) (QR, start mikrokanoniczny, A′) | p-Ps 137,4 ps, o-Ps 150,1 ns (audyt 359) | 125,14 ps / 142,04 ns | \(+9{,}8\%\) / \(+5{,}6\%\); rozpad niejednowykładniczy (rozrzut \(L\) ze startu) |
| \(\tau_{\rm orto}/\tau_{\rm para}\) | \(1113{,}97\) (start Langera); 1114 (mikrokanoniczny, ten sam skład prób); 933–953 bez transportu singletowego | \(1135{,}0\) | reguła wyboru \(\varepsilon_{\rm OP}\) jest importem; własny wkład modelu to \(w\) |
| \(\tau_{\rm orto}/\tau_{\rm para}\), eksperyment 6 (stan końcowy, audyt 375) | silnik \(1/\varepsilon=1113{,}9\) (124,494 ps / 138,674 ns); próba 1000 par \(1248\pm4{,}5\%\), rozkłady wykładnicze (KS \(p=0{,}90\) / \(0{,}99\)) | \(1135{,}0\) | \(-1{,}9\%\) (\(\varepsilon\) to QED wiodącego rzędu); stary mechanizm przy barierze (`--contact-barrier`) dawał \(27\,215\pm42\%\), z klasycznym czynnikiem geometrii \(24{,}4\) |
| kwantyzacja spinów przez pole ziemskie (50 µT, spiny swobodne, stan końcowy \(n=1\), \(L=\hbar/2\); audyt 376) | \(\tau(50\,\mu{\rm T})/\tau(0)=1{,}0000002\) w obu kanałach; mediana \(\max_t|\Delta w|\) do 1 µs \(8\cdot10^{-6}\) (p-Ps) / \(1\cdot10^{-5}\) (o-Ps); liczba par z \(w<0{,}01\) lub \(>0{,}99\) bez zmian (10 / 3 z 200); o-Ps zostaje przy \(1/\langle\Gamma\rangle=0{,}35\) ns | o-Ps 142,04 ns | pole **nie kwantuje**: precesja przeciwbieżna e⁻/e⁺ (\(8{,}8\cdot10^6\) rad/s) jest uśredniana przez wspólny obrót \(\mathbf L\cdot\mathbf S\) wokół \(\mathbf J\) (\(\sim1{,}65\cdot10^{12}\) rad/s); kontrola 1 T: \(w\) miesza się (chaotycznie), kanały zbliżają się (p-Ps \(\tau\times1{,}12\), o-Ps \(\times0{,}96\) — znak jak w gaszeniu magnetycznym), ale bez przyciągania do \(w=0/1\); rzut \(\mathbf J\) na pole zachowany do \(10^{-11}\hbar\) z ciągłą wartością początkową, a rzut momentu też ciągły — brak kwantyzacji przestrzennej (376f). Kwantowo moment Ps jest zerowy we wszystkich stanach (brak liniowego Zeemana); w modelu o-Ps ma 0, a klasyczny p-Ps \(2\mu_B\) |
| linia 1S–2S do rzędu \(\alpha^4\) (różnica poziomów \(1^3S_1\), \(2^3S_1\)) | \(+96{,}8\) (relatywistyczna) \(-25{,}6\) (kontakt) \(-51{,}2\) (pętla) GHz względem \(\tfrac34B\): 5,102218 eV (audyt 365) | 5,101790 eV (Fee i in. 1993) | \(+103{,}5\) GHz (bez pętli \(+154{,}7\)); brak klasycznego odpowiednika wirtualnej anihilacji (kwantowo \(-76{,}6\) GHz) |
| czas życia 2p / 3p wodoru (`--pair proton,electron`, start Langera) | **1,681 / 5,675 ns** (domyślnie, audyt 370: tempo \(=\dot L_{\rm klas}/\hbar\); z tablicą \(S(e)\) 1,665 / 5,428 ns; z hazardem z przerwy 4,995 / 10,18 ns). Moc silnika = Larmor do \(10^{-4}\); +5,1% to dokładność przepisu półklasycznego: z \(L=\sqrt{l(l+1)}\hbar\) wychodzi \(-6{,}6\%\) / \(-9{,}6\%\) (audyt 371) |
| czas życia 4p / 5p wodoru | Langer 13,45 / 26,27 ns (\(+19{,}6\%\) / \(+20{,}0\%\)); \(\sqrt{l(l+1)}\): 11,96 / 23,36 ns (audyt 372) | 11,25(78) / 21,9(3,0) ns (Etherton i in., Phys. Rev. A 2, 2177 (1970)) | błąd **nie maleje** z \(n\) przy stałym \(l=1\) (orbita coraz bardziej radialna); żadne jedno przypisanie \(L\) nie pasuje do całego ciągu | 1,600(4) / 5,58(13) ns (Bickel i Goodman, Phys. Rev. 148, 1 (1966)) | klasyczne tempo \(S(e)P_{E1}/\hbar\omega_{\rm orb}\) zgodne do \(+4{,}1\%\) / \(-2{,}7\%\); hazard z przerwy (od 352) spowalnia o \(E_{\rm przerwa}/\hbar\omega\) (3 / 1,875); pomiaru 2P pozytonium brak |
| drgania w polu punktu zerowego (Welton, widmo SED modelu) | stany S w górę: 1S 5,72, 2S 0,86 GHz przy paśmie \([\omega_n, mc^2/\hbar]\); 1,25 / 0,16 GHz przy paśmie modelu \([0{,}3;3]\omega\) (audyt 373) | — | właściwy znak i skala \(\alpha^5\ln\), ale tylko ~5% pozostałych niezgodności (linia 1S–2S 103,5 → 98,7 GHz; \(2^3S_1\) 14,8 → 13,95 GHz); mody ZPF modelu zaniżają wariancję drgań (0,59–0,94) |
| nakładające się pętle przy kontakcie (dokładnie) | magnetycznie \(1{,}023\pm0{,}009\times\) kontakt Fermiego, elektrycznie \(1{,}030\pm0{,}015\times\) człon skończonego rozmiaru (audyt 366) | — | dokładnie oba człony, które model już ma; **brak członu tylko trypletowego** — wirtualna anihilacja (zniknięcie i powstanie pary) nie ma klasycznego odpowiednika w modelu z zachowaną liczbą cząstek |
| skończony rozmiar ładunku (pętla \(R=2r^*\)) | stany S \(+\alpha^4mc^2(g/2)^2/(6n^3)\): 1S 58,53, 2S 7,32 GHz (audyt 365) | — | bez nowego parametru; \(1{,}34\times\) kwantowy człon Darwina; deficyt \(2^3S_1\) 22,12 → 14,80 GHz |
| odstęp \(2^3P_1-2^3P_0\) | 2,98 GHz (V1); 3,04 / 2,84 GHz w innych klasycznych modelach wektorowych (audyt 374) | 5,487 GHz | \(-45\) do \(-48\%\) w każdym przypisaniu: klasyczna energia tensorowa J = 0 (S antyrównoległe do L) nie różni się od J = 2; kwantowe \(\langle S_{12}\rangle=-4\) dla \(^3P_0\) wymaga koherentnej superpozycji składowych spinu — ograniczenie strukturalne; \(P_2-P_1\) zgodne do 1% |
| \(2^3S_1-2^1P_1\) | \(-10{,}94\) GHz (audyt 364) | \(+11{,}180\) GHz (Ley i in., Hyperfine Interact. 89, 327 (1994)) | izoluje stan S: \(2^3S_1\) za nisko o 22,12 GHz; część spinowa \(2^3P_J\) myli się tylko o \(-0{,}7\) GHz (J = 1, 2) i \(-3{,}3\) GHz (J = 0) |
| struktura subtelna \(n=2\) (\(2^3S_1-2^3P_{0,1,2}\)) | \(-6{,}88\), \(-9{,}86\), \(-14{,}19\) GHz (audyt 363) | \(+18{,}500\), \(+13{,}012\), \(+8{,}624\) GHz (Hagena i in., PRL 2887 (1993)) | zły znak: \(2^3S_1\) poniżej \(2^3P_J\); brakuje wirtualnej anihilacji i kontaktowego członu Darwina w stanach S (16,4 GHz) oraz ~6–9 GHz części spinowej P; \(P_2-P_1\) zgodne do 1,3% |
| rozszczepienie nadsubtelne 1s | \(58{,}53\) GHz (\(n=1\), \(L=\hbar/2\), kontakt QR, audyt 360); dawniej: koło \(L=\hbar\) \(6{,}7\cdot10^{-7}\) pomiaru, zespół mikrokanoniczny \(290{,}5\) GHz (324) | \(203{,}3942\) GHz | \(0{,}288\) pomiaru, znak poprawny; \(=\tfrac12\cdot\tfrac47\) wartości wiodącego rzędu: \(\tfrac12\) z klasycznych korelacji spinów, \(\tfrac47\) bo wirtualnej anihilacji (3/7) klasycznie nie ma |
| rozpraszanie | kształt Rutherforda, jedna normalizacja | — | klasyczne |

**Wynik negatywny jest wynikiem.** Różnica para/orto w przyrodzie pochodzi z
reguły wyboru \(2\gamma/3\gamma\) (parzystość ładunkowa), której klasyczna
dynamika nie zawiera; w modelu wchodzi przez \(\varepsilon_{\rm OP}\) i
\(w\). Skala bezwzględna wymaga gęstości w zerze, której klasyczna orbita
nie ma; tożsamość QR daje ją z \(\langle dV/dr\rangle\), ale jedno \(L\) na
parę daje rozrzut temp.

## 10. Granice dziedziny i otwarte problemy

- **Bariera Comptona** \(r^*\): poniżej niej klasyczna elektrodynamika
  punktowa przestaje obowiązywać, a przebieg się kończy.
- **Granica retardacji:** okres / czas przelotu światła \(\le150\) — margines
  przybliżenia Darwina i całkowania, nie fizyka.
- **Poniżej \(n=1\):** w trybie fotonowym z regułą działania nie ma
  poziomów; w trybie ciągłym orbita przechodzi klasyczną katastrofę
  promienistą (\(\approx199\) ps), której prawdziwe pozytonium nie przechodzi.
- **Otwarte** (lista z audytów do 359, uaktualniona w audycie 385):
  (1) jedno \(L=\hbar/2\) na \(n=1\) daje tylko przygotowanie Langera
  (import WKB, nr 6) — dynamika go nie wybiera; audyt 386: to wybór ścisłego
  \(\langle r^{-2}\rangle\), a żaden rozkład \(L\) przy stałym \(a\) nie odtwarza
  stanu 1s, więc test dynamiczny (SED) dotyczy rozkładu \((E,L)\); audyt 387: w SED
  (Fokker–Planck w zmiennych działania) orbita \(n=1\) nie ma stanu
  stacjonarnego — 90% się jonizuje, 10% kolapsuje — więc SED niczego nie
  wybiera; zamknięte negatywnie; (2) oscylacja \(w\) w p-Ps
  jest usunięta transportem singletowym (import nr 8, audyt 355), a dryf
  \(L\) w o-Ps (356) był błędem całkowania, usuniętym w audycie 358.
  Dawne punkty zostały rozstrzygnięte: tempo kroku \(1\to0\) jest w silniku
  od audytu 361 (sekcja 6.2); przy \(|L|<\hbar\) obowiązuje \(\Delta l=+1\)
  (audyt 361, import nr 3); czasu kaskady \(2\to1\) nie da się porównać z
  pomiarem \(2P\) pozytonium, bo takiego pomiaru w repozytorium nie ma (dla
  wodoru: audyty 370–372).
- **Nie ma w modelu:** funkcji falowej, zasady Pauliego, poprawek
  radiacyjnych, rozpraszania Bhabhy, odpowiedzi detektora ani zjawisk w
  materiale (pełna lista w README, „Ograniczenia”).

- ~~Kontakt Fermiego w transporcie spinów i w bilansie energii spinów nie był
  gęstością QR~~ — **rozwiązane w audycie 378.** Audyt 377 zmierzył kontakt
  Plummera z orbity: \(6{,}4\)–\(7{,}4\cdot10^{-4}n_{\rm QR}\) przy \(l=0\),
  a przy \(l\ge1\) niezerowy. Teraz model ma jedną gęstość (sekcja 3).
  Czasy życia się nie zmieniły, a linia \(2p\to1s\) rozdzieliła się o HFS 1s.

- **Odrzucanie fotonów przy \(k\ge2\)** (audyty 381, 383). Silnik daje
  \(\Delta n=\min(k,n-1)\) i odrzuca foton, gdy \(L'\) nie mieści się w
  \(n'\). Nie skraca skoku do dozwolonego. Czasy życia wodoru (rachunek
  półanalityczny, audyt 383) wobec QM (teoria, nie pomiar):

  | stan | obecnie (odrzucanie) | skracanie skoku | QM |
  |---|---|---|---|
  | 3d / 4d / 4f | 31,8 / 129 / 131 ns | 15,77 / 37,4 / 73,2 ns | 15,6 / 36,5 / 72,6 ns |
  | 2p–5p | bez odrzuceń | 1,68 / 5,68 / 13,45 / 26,3 ns | 1,60 / 5,40 / 12,4 / 23,6 ns |
  | 3s / 4s | 145 / 120 ns (tablica), 104 / 183 ns (dokładne \(w_k\)) | 0,63 / 1,5 ns | 158 / 226 ns |

  - **Stany z \(l\ge1\)** (\(\Delta l=-1\)): skracanie skoku zgadza się z QM
    do 1–8%, a odrzucanie zawyża czasy stanów d i f 2–3,5 razy. Stany p się
    nie zmieniają: to one mają pomiar (2p–5p), a wszystkie ich fotony są
    przyjmowane.
  - **Stany s** (\(\Delta l=+1\)): klasyczny moment siły promieniowania
    zawsze zmniejsza \(L\), więc foton zwiększający \(L\) nie ma klasycznego
    źródła, a hazard \(|dL/dt|/\hbar\) liczy fotony niedostępne dla stanu s.
    Udział harmonicznej \(k=1\) (zasada korespondencji) jest zastępczym
    oszacowaniem: \(-34\%\) / \(-19\%\) wobec QM przy dokładnych \(w_k\).
    Amplituda Kramersa z orbity pośredniej daje 30,6 / 60,8 ns, czyli jest
    gorsza.
  - **Tablica udziałów harmonicznych** (`eccentricOrbitHarmonicNumber`,
    \(e\le0{,}98\), interpolacja kwantyli) myli \(w_1\) o \(+26\%\) przy
    \(e=0{,}866\), \(-38\%\) przy 0,968 i \(-30\%\) przy 0,986 (dokładnie:
    \(w_k\propto k[J_k'(ke)^2+\tfrac{1-e^2}{e^2}J_k(ke)^2]\)).
  - **Wprowadzone w audycie 384:** skracanie skoku dla \(\Delta l=-1\)
    (`CREM_REFUSE_ACTION_STEP=1` przywraca odrzucanie) i dokładne \(w_k\)
    z funkcji Bessela (`CREM_HARMONIC_TABLE=1` przywraca tablicę). Silnik
    daje dla wodoru 3d pierwszy foton po 15,764 ns (wcześniej ~31,8; QM
    15,6 — teoria), 2p i 3p bez zmian (1,681 / 5,675 ns). **Otwarte
    zostają stany s:** foton \(\Delta l=+1\) bez klasycznego źródła, oceniany
    przez udział \(k=1\). Dla wodoru 3s / 4s daje to 104 / 183 ns wobec
    teorii QM 158 / 226 ns.
- **Koszt transportu spinów w stanach s** (audyty 381, 385). W 3s Ps
  transport spin–orbita potrzebuje ~773 podkroków na pół checkpointu
  (1,26 ns symulacji ≈ 96 s zegara). Próba z audytu 385 — dokładny obrót
  \(R_J(A|J|t)\,R_S((\kappa-A)|S|t)\) dla części \(\mathbf L\cdot\mathbf S\)
  i Heisenberga — zachowuje kąt między spinami do \(10^{-11}\), ale reszta
  nakładana połówkami kroku wokół obrotu do \(\pi\) psuje orientację
  (2,9e-4 w jednym kroku, O(1) po 100 ps). Zbiega dopiero przy kroku ~0,05
  rad, a wtedy zysk znika. Nie jest domyślna (`CREM_SPIN_S_STATE_EXACT=1`).
  **Otwarte:** reszta w obrazie oddziaływania.
  Audyt 388: jedna średnia BMT zamiast dwóch daje 1,23× przy identycznych
  wynikach; ~80% kosztu podkroku przy 3s to magnetyczna precesja
  perycentrum (siły z momentami w 512 węzłach).
- ~~Zegar laboratoryjny w długich przebiegach~~ — **naprawione w audycie
  382.** Czas laboratoryjny był mnożony przez \(\gamma\) dryfu dipolowego
  (audyt 110). Po audytach 124–127 przyrost pędu na orbitę to reszta
  całkowania, ale silnik sumował go liniowo przez \(\sim10^5\) orbit na
  checkpoint: w 3s \(\beta\approx0{,}7\) po 6 ns, a potem \(10^{153}\) ps.
  Zegar laboratoryjny uwzględnia teraz tylko odrzut fotonów, a dryf zostaje
  diagnostyką w wydruku (`CREM_LAB_DIPOLE_DRIFT=1` przywraca audyt 110).
  Przy \(n=1\) różnica wynosiła \(\gamma-1=4\cdot10^{-9}\), poniżej
  drukowanej dokładności.
- **Koszt stanów s powyżej \(n=1\)** to ~70–100 s zegara na 1 ns symulacji
  (spiny: ~773 podkroki na pół checkpointu). Przy ~207 ns życia 3s Ps (od audytu 384; wcześniej
  ~600 ns) daje to kilka godzin na trajektorię.

- **Wychwyt e⁺e⁻ jest klasyczny** (audyty 393–394). Eksperyment 5 daje
  przekrój wychwytu zgodny z klasycznym wychwytem radiacyjnym (energia
  wypromieniowana przez dipol na jeden przelot > E): 2,43·10⁷ b wobec
  2,82·10⁷ b przy 1,25 eV, ~3700 razy więcej niż Kramers, który liczy
  wychwyt jednym fotonem o energii > E. **Rozwiązane w audycie 395:** w
  trybie `stochastic` interferencyjna część promieniowania pary była liczona
  dwa razy (w fotonach z pełnej mocy dipola i w ciągłej wymianie przez
  wzajemne pola opóźnione). Teraz siła −M_i = −(q_i/q_j)·LL_j kasuje wiodący
  człon tej wymiany (ciągła strata 0,56 → 0,06 mocy Larmora pary), a wychwyt
  w eksperymencie 5 przy `--radiation-reaction stochastic` spada z 19 do 0
  na 400 zdarzeń (Kramers ~0,1). Czasy życia Ps bez zmian.

- **Spłata fotonu na ścieżce mechanicznej** (audyty 398–399). Foton pary
  związanej niesie tam E(n) − E(n − Δn) (reguła działania, 398). Od audytu 399
  okno spłaca też ΔJ = Δn ħ i ΔL = mħ (reguła Δl: E1 axialPhotonSpinSign, E2
  m = +2 przez photonOrbitalTransferHbar). J i L idą w jednym tempie, energia
  przyrostu to E(J) − E(J − dJ), a L ustawia obrót pędu względnego w układzie
  CM przy stałym |p| (nic innego się nie zmienia) — do celu bezwzględnego
  L_kon + ρħ·(zaległe J), ρ = ΣΔL/ΣΔJ (audyt 400), więc L odebrane ponad
  miarę przy apsydzie (L ≤ r|p⊥|) wraca tam, gdzie p_r ≠ 0. Cel nieosiągalny w
  danym miejscu (|L_cel| > r|p⊥|) zostawia pęd bez zmian (audyt 401): przycięcie
  do r|p⊥| zerowało p_r, robiło z każdego punktu apsydę i przyszpilało r
  (`CREM_LADDER_CLAMP_ROTATION=1` przywraca). Okno jest adiabatyczne:
  trwa max(1, ΔJ/(ε n)) okresów, ε = 0,01, bo okno jednego okresu wysysało
  energię kinetyczną przy stałym r i L spadało do 0,04ħ. Foton odpalony przy
  otwartym oknie liczy lukę od n − (zaległe J). Sonda: z kołowej n = 2 okno
  kończy na n = 0,99993, L = 0,9988ħ, e = 0,047 i orbita stoi (398: zapaść do
  n = 0,09); od audytu 400: n = 0,999986, L = 0,99936ħ, e = 0,035. Dwa
  nakładające się okna (3 → 2 → 1, artefakt hazardu ×3e5) kończą na
  n = 1,00008, L = 1,000003ħ, e = 0,012 (w 399, z celem przyrostowym: n =
  1,007, e = 0,19 — L przeciągnięte o 1,16ħ przy apsydach, orbita e ≈ 0,77, na
  której dynamika dodała 0,093 eV; teraz 1,7e-22 J). W połowie złożonego okna
  e chwilowo rośnie (L = 1,22ħ przy n = 1,97) i wraca; od audytu 401 max e w
  oknie 0,032, koniec n = 1,000026, L = 1,000018ħ. Sprawdzone też (401): orbita
  e₀ = 0,6 (okno zachowuje J_r = n − L/ħ: 0,600 → 0,604) i E2 w wodorze
  (wymuszone, n = 3 → 1, foton 12,09 eV, koniec L = 0,9995ħ, e = 0,029).
  `CREM_WINDOW_SCALE_L=1` przywraca skalowanie (398).
- **Dryf energii ścieżki mechanicznej na ciasnych orbitach ekscentrycznych**
  (audyt 401). Bez fotonów, n = 1, 300 obiegów, tryb `stochastic` (chwilowe
  siły Coulomb–Darwin + dipole między fotonami): e = 0 → Δn = +0,0019/1500
  obiegów, e = 0,5 → +0,068/1500, e = 0,917 → +0,047/300. To błąd całkowania:
  w energii zachowanej taki sam, rzędu 1 w kroku zewnętrznym (128/512/1024/2048
  kroków na obieg: +0,0467/+0,0153/+0,0072/+0,0036), zależny od tolerancji
  (1e-7: +0,0082), ~55 % bez dipoli, niezależny od prawa kroku (stały kąt:
  +0,0487). W trybie z polami opóźnionymi ten sam
  zysk maskuje fizyczną stratę (`disabled`: +0,0038 przy 128, −0,0236 przy 512
  krokach). **Przyczyna (audyt 402):** kick–drift–kick liczy siły końca kroku
  w stanie z prędkością po POŁOWIE kopnięcia; człony zależne od prędkości
  (Darwin, v × B dipoli, ładunek–dipol) dostają błąd lokalny O(dt²), globalny
  rzędu 1. Jedna ponowna ewaluacja przy prowizorycznej prędkości końcowej
  (`CREM_VELOCITY_CORRECTOR=1`, test) daje: krok 128 +0,0178, 512 +0,00017
  (rząd ~3,4), udział dipoli znika, `disabled` −0,0251 już przy 128. Wpływ na
  eksp. 5 (400, ziarno 394, tryb domyślny): związane 41 → 38, σ_capture
  3,00 → 2,76·10⁷ barn (±0,45) — w szumie; czas ×1,09. W trybie `stochastic`
  korektor dodaje 3 błędy numeryczne przy przelotach blisko bariery (r_min
  0,2–0,4 pm). Audyt 403: to odrzucenia dokładności na maksymalnej głębokości
  tuż nad progiem (1,05–1,56e-5 przy 1e-5) w zdarzeniach, które i bez
  korektora były marginalne; eksp. 5 całkuje z polami opóźnionymi, a w nich
  siła gradientowa dipola (covariantDipoleGradientForce) skacze o 60 % między
  połowieniami dt (nieciągłość C0, audyt 149), sektor Coulomba jest czysto
  drugiego rzędu. **Przyczyna (audyt 404):** moment siły reakcji M1,
  m × (m₁‴ + m₂‴)·μ₀/(6πc³) z m‴ różnicowanym jednostronnie po historii
  (spinowy odpowiednik Abrahama–Lorentza całkowany jawnie), wpada przy
  barierze w sprzężenie zwrotne: 3,8·10²¹ 1/s wobec precesji BMT 2,4·10¹⁹
  i fizycznego ~10¹⁴; rozrzuca momenty o ~20 % między węzłami historii, m″
  źródła wychodzi 100–2000 × ω²m, pola dipola są zdominowane przez ten
  fałszywy człon promieniowania (E ≈ 15 × Coulomb), a siła gradientowa jest
  resztą 10⁻⁹ dwóch członów ~10¹⁶. `CREM_M1_TORQUE_REDUCED=1` (test) bierze
  m‴ = Ω × (Ω × (Ω × m)) z precesji BMT (jak Landau–Lifszyc dla ładunku):
  1,4·10¹⁴ 1/s; eksp. 5 (400, ziarno 394) bez cenzury i błędów (domyślny
  400/0/0 zamiast 388/12/0, szybciej; stochastic z korektorem 400/0/0 zamiast
  387/10/3, 2 wychwyty Ortho z dawniej ocenzurowanych przelotów przy
  barierze), walidacja 65/65, Ps bez zmian. **Otwarte:** promocja obu
  (zredukowany moment M1 i korektor) do domyślnych; moc M1 nadal z m″ historii.
  Pomiar: `CREM_MECHANICAL_STEPS_PER_ORBIT`, `CREM_MECHANICAL_TOLERANCE`.

## 11. Gdzie szukać

| temat | miejsce |
|---|---|
| prawa sił, pola, rozmycie dipola | `modules/electrodynamics.hpp` |
| integrator, historia retardowana | `modules/crem_engine.hpp` |
| przygotowanie stanu, pętla trajektorii | `modules/crem_trajectory.hpp` |
| tryb fotonowy, reguła działania, sufit, zamknięcie, tempo anihilacji | `modules/crem_collapse.hpp` |
| transport \(\mathbf L,\mathbf S_1,\mathbf S_2\) | `modules/secular_spin_orbit.hpp` |
| eksperyment 6, fotony anihilacji | `modules/contact_annihilation.hpp` |
| wszystkie przełączniki i ich uzasadnienia | `modules/configuration_panel.hpp` |
| pomiary, wyceny importów, historia | `README.md` |
| rozumowanie krok po kroku, korekty | plik audytu i `audits/data/` |
