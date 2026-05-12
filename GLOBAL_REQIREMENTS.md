Пишем набор скриптов. Задача: обработать огромное количестве (1.2 млн строк) данных о трёхмерных белковых струкутрах из PDB, расширив описание.

Это - отдельные структурные мотивы и/или домены, найденные в белках. Наша задача - расширить описание этих структур, пересчитать ряд геометрических параметров и охарактеризовать их.

# Предусловия:

- файлы PDB лежат тут: `/home/gluck/mnt/Data2/PDB/pdb/`
- файлы отсортированы по папкам с двумя средними символами кода PDB, пример пути к структуре: `/home/gluck/mnt/Data2/PDB/pdb/01/pdb201d.ent.gz`
- Они, как можно увидеть, в .gz-архивах
- Файл данных лежит тут: `data/data-1778568752028.csv`
- есть файл _HEAD с первой тысячей строк: `data/data-1778568752028_HEAD.csv`
- Реально работаем с файлом `data/data_SORTED.csv`, от отсортирован по experiment_id, chain, motif_start

Результат `print(df.describe().transpose().to_markdown())`

|                |   count |             mean |            std |            min |             25% |             50% |              75% |              max |
|:---------------|--------:|-----------------:|---------------:|---------------:|----------------:|----------------:|-----------------:|-----------------:|
| id             |     999 |      4.47013e+06 | 528453         |    3.86794e+06 |     4.15996e+06 |     4.41653e+06 |      4.69158e+06 |      1.08004e+07 |
| model_id       |     999 |      0.994995    |      0.0706041 |    0           |     1           |     1           |      1           |      1           |
| chain_len      |     999 |      0.706707    |     10.0986    |    0           |     0           |     0           |      0           |    172           |
| resolution     |     999 |      3.04883     |      4.37598   |    0           |     2           |     2.5         |      3.1         |     70           |
| motif_len      |     999 |     42.02        |     13.4624    |   16           |    33           |    39           |     50           |    115           |
| motif_start    |     999 |    306.715       |    371.805     |    0           |    86           |   203           |    377           |   3891           |
| motif_end      |     999 |    350.034       |    372.617     |   21           |   132.5         |   245           |    423.5         |   3929           |
| organism_id    |     999 | 152415           | 301021         |    0           |  9606           | 10090           | 208964           |      2.82942e+06 |
| rg             |     999 |     11.7178      |      2.72678   |    6.84892     |     9.91561     |    11.1088      |     13.1162      |     29.9358      |
| sasa           |     999 |   3754.54        |   1056.49      | 1698.74        |  2994.53        |  3609.14        |   4342.64        |   9177.6         |
| surface_k      |     999 |      0.546094    |      0.149304  |    0.104326    |     0.446448    |     0.552684    |      0.643252    |      1           |
| inner_hb       |     999 |      9.05806     |      5.14155   |    0           |     5           |     8           |     12           |     37           |
| inner_hp       |     999 |     33.8969      |     21.0015    |    2           |    21           |    29           |     42           |    270           |
| inner_ds       |     999 |      0.011011    |      0.104406  |    0           |     0           |     0           |      0           |      1           |
| outer_hb       |     999 |      8.06306     |      4.87037   |    0           |     5           |     7           |     11           |     28           |
| outer_hp       |     999 |     45.9249      |     26.228     |    0           |    30           |    42           |     59           |    306           |
| outer_ds       |     999 |      0.043043    |      0.267005  |    0           |     0           |     0           |      0           |      3           |
| motif_weight   |     999 |   4718.29        |   1530.75      | 1905.1         |  3584.29        |  4458.54        |   5641.34        |  12497.4         |
| protein_weight |     999 |  74282.9         |  73123.2       | 6025.61        | 33326.6         | 53831.9         |  89878.3         | 837769           |
| weight_percent |     999 |     11.3286      |     10.4836    |    0.363743    |     4.60121     |     7.99723     |     13.4572      |     98.0576      |
| binding_sites  |       0 |    nan           |    nan         |  nan           |   nan           |   nan           |    nan           |    nan           |
| helix_count    |     999 |      1.74575     |      0.652895  |    0           |     2           |     2           |      2           |      4           |
| strand_count   |     999 |      0.525526    |      1.29787   |    0           |     0           |     0           |      0           |      7           |
| coil_count     |     999 |      3.25225     |      0.762624  |    1           |     3           |     3           |      3           |     10           |
| helix_length   |     999 |     27.0811      |     16.0394    |    0           |    18           |    27           |     35           |    107           |
| strand_length  |     999 |      2.95495     |      7.70218   |    0           |     0           |     0           |      0           |     44           |
| coil_length    |     999 |     11.984       |      7.2337    |    3           |     8           |     9           |     13           |     55           |

Файл является выгрузкой из БД (PostgreSQL) следующей структуры:

```sql
-- Table: public.motif

-- DROP TABLE IF EXISTS public.motif;

CREATE TABLE IF NOT EXISTS public.motif
(
    id integer NOT NULL DEFAULT nextval('motif_id_seq'::regclass),
    experiment_id character varying(32) COLLATE pg_catalog."default" NOT NULL,
    source character varying(16) COLLATE pg_catalog."default" NOT NULL,
    uniprot_id character varying(16) COLLATE pg_catalog."default" NOT NULL,
    protein_name character varying(256) COLLATE pg_catalog."default",
    model_id integer NOT NULL,
    chain_id character varying(2) COLLATE pg_catalog."default" NOT NULL,
    chain_len integer NOT NULL,
    organism character varying(256) COLLATE pg_catalog."default",
    resolution real,
    experiment_type character varying(16) COLLATE pg_catalog."default",
    experiment_name text COLLATE pg_catalog."default",
    biological_process character varying(1024) COLLATE pg_catalog."default",
    molecular_function character varying(1024) COLLATE pg_catalog."default",
    fasta text COLLATE pg_catalog."default" NOT NULL,
    apfid character varying(64) COLLATE pg_catalog."default" NOT NULL,
    motif_type character varying(32) COLLATE pg_catalog."default" NOT NULL,
    motif_len integer NOT NULL,
    motif_start integer NOT NULL,
    motif_end integer NOT NULL,
    full_chain boolean DEFAULT false,
    single_chain boolean DEFAULT false,
    organism_id character varying(16) COLLATE pg_catalog."default",
    rg real,
    sasa real,
    surface_k real,
    inner_hb integer,
    inner_hp integer,
    inner_ds integer,
    outer_hb integer,
    outer_hp integer,
    outer_ds integer,
    motif_weight real,
    protein_weight real,
    weight_percent real,
    protein_partners text COLLATE pg_catalog."default",
    has_binding_sites character varying(3) COLLATE pg_catalog."default",
    binding_sites text COLLATE pg_catalog."default",
    helix_count integer,
    strand_count integer,
    coil_count integer,
    helix_length integer,
    strand_length integer,
    coil_length integer,
    CONSTRAINT motif_pkey PRIMARY KEY (id)
)

TABLESPACE pg_default;

ALTER TABLE IF EXISTS public.motif
    OWNER to data_api;
-- Index: motif_apfid

-- DROP INDEX IF EXISTS public.motif_apfid;

CREATE UNIQUE INDEX IF NOT EXISTS motif_apfid
    ON public.motif USING btree
    (apfid COLLATE pg_catalog."default" ASC NULLS LAST)
    TABLESPACE pg_default;
-- Index: motif_experiment_id

-- DROP INDEX IF EXISTS public.motif_experiment_id;

CREATE INDEX IF NOT EXISTS motif_experiment_id
    ON public.motif USING btree
    (experiment_id COLLATE pg_catalog."default" ASC NULLS LAST)
    TABLESPACE pg_default;
-- Index: motif_resolution

-- DROP INDEX IF EXISTS public.motif_resolution;

CREATE INDEX IF NOT EXISTS motif_resolution
    ON public.motif USING btree
    (resolution ASC NULLS LAST)
    TABLESPACE pg_default;
-- Index: motif_uniprot_id

-- DROP INDEX IF EXISTS public.motif_uniprot_id;

CREATE INDEX IF NOT EXISTS motif_uniprot_id
    ON public.motif USING btree
    (uniprot_id COLLATE pg_catalog."default" ASC NULLS LAST)
    TABLESPACE pg_default;
```

К БД НЕ обращаемся и не пытаемся, это просто для примера. Выходной формат - такой же, как входной.

# Требуемые данные

Первое, самое важное: apfid

id - просто ID в базе данных. В принципе игнорируем.

apfid - это идентификатор. Подробнее можно прочитать тут: https://raw.githubusercontent.com/protdb/apfid/refs/heads/master/README.MD

Он остаётся у всех структур, по нему же данные будут обновляться в результате.

Параметры структуры:
- model_id: ID модели в PDB-файле ,если их несколько. Считаем, что во всех случаях используется самая первая модель. В данных встречается то 1, то 0, пока игнорируем это различие и оставляем as is.
- chain_id - ID цепи (буква)
- chain_len - длина цепи в аминокислотных остатках. Необходимо пересчитать
- motif_len - длина самого мотива (в аминокислотных остатках)
- motif_start, motif_end - начало и конец
- rg: радиус гирации/радиус вращения. Нужно рассчитать
- sasa: SASA, необходимо рассчитать. Величина SASA мотива, взятого отдельно от остальной белковой цепи.
- surface_k - Коэффициент погруженности мотива. Для этого параметра считаем SASA для мотива в белке (т.е. SASA только аминокислот мотива) и делим на SASA мотива вне белка, т.е. без учета влияния других аминокислот белка
- поля inner/outer hp/hb/ds - взаимодействия (hp - гидрофобные, hb - водородные связи, ds - дисульфидные взаимодействия) внутри (inner) мотива (когда обе взаимодействующие аминокислоты находятся внутри мотива) и внешние (outer), когда взаимодействие между аминокислотой мотива и внешней (в цепи)
- motif_weight - молекулярная масса мотива
- helix/strand/coil count - количество элементов вторичной структуры в мотиве, рассчитывается с использованием STRIDE, учитываются законченные структурные элементы в рамках мотива (скажем, если мотив состоит из двух спиралей и перетяжки между ними, то helix_count=2, coil_count=1; которкие coil'ы на концах структуры после обрезки не учитываются)
- helix/strand/coil length - суммарная дилна элементов вторичной струкутры

Добавляем параметр: алифатический индекс (GRAVY). Расчет на основе алгоритма здесь: `/home/gluck/PycharmProjects/CornerPhobia/corner_hydropathy/`. Нужно учесть, что нам отсюда нужен только один параметр, остальные слишком специфические.


Пареметры белка получаем через uniprot_meta_tool, в соответствии с примером здесь: `examples/get_up_additionals.py`. Документация по библиотеке и исходный код есть тут: `/home/gluck/PycharmProjects/uniprot_infotool/README.md`

Используем поле uniprot_id из исходных данных.

Получаем:
- protein_name: перезаписываем
- biological_process
- molecular_function
- protein_weight: пересчитываем из fasta белка
- weight_percent: motif_weight / protein_weight * 100
- protein_partners
- has_binding_sites
- binding_sites

Большая часть есть в примере либо в исходниках библиотеки

# Расчет данных

Не строим лишние велосипеды!

- Используем biopython для основной массы операций
- Используем FreeSASA для расчетов SASA
- Для расчета взаимодействий используем ProDy
- Обработку STRIDE строим на основе скрипта `examples/stride.py`, но подходим к нему критически, делаем новый и переписываем по необходимости
- Молекулярный вес рассчитываем с biopython. Для последовательностей с X/UNK считаем среднюю массу аминокислотного остатка

# Проверки данных:

Могут быть проблемы с тем, что FASTA струкутры в БД не соответствует реальности. Могут быть координаты некорректные координаты, по которым вырезается пустая структура. Такие струкутры нужно сохранять в отдельный список на удаление из БД.

Нужен отдельный флаг конфигурации STOP_ON_ERROR (см. ниже про конфигурацию). Если False - то пишем ошибки в лог в формате:
```
==========================================================
{apfid}
{traceback}
```

# Конфигурация скрипта

- Выносим в отдельный файл config.py констатнами все пути к файлам, STOP_ON_ERROR, а также START_INDEX (в файле, откуда начинать обработку) и STOP_INDEX - завершать на нём или смотреть дальше

# Описание процесса

В общих чертах так:
- Читаем файл PDB
- Выполняем расчеты для файла
- Вырезаем мотивы по координатам (строго первая модель, игнорируем остальные)
- Сохраняем в файл, определяем вторичную струкутру через STRIDE
- Выполняем остальные вычисления
- Получаем данные UniProt, получаем данные по белку и мотиву
- Сохраняем во временный файл

# Технические требования

Процесс долгий и ресурсоёмкий, также возможны остановки и перебои в выполнении. Поэтому:
- Выполняем расчеты многопоточно, используем multiprocessing. обработка одной струкутры - линейная, выполняем её одним процессом
- Результаты по каждой отдельной структуре сохраняются в отдельную папку, результаты расчетов сохраняем в отдельную папку, по одной струкутре в json-файл. в имени файла сохраняем apfid, таким образом список обработанных файлов можно получить через список файлов в папке
  - также сортируем по подпапкам с именами 0000, 1000, 2000 и т.д. - по порядку следования структуры в исходной таблице. сопоставляем только ID файлов в текущей задаче (START_INDEX -> STOP_INDEX), чтобы не листать всю папку на этапе, когда готовы десятки тысяч струкутр.
- То же касается некорректных структур (пустые координаты): складываем в ту же папку, но другая струкутра JSON, в имени после apfid добавляем __incorrect
- Добавляем флаг ONLY_NOT_PROCESSED в конфигурации: если True, то исключаем отмеченные
- Место на диске экономим. вырезаем PDB-струкутру, обрабатываем через STRIDE, сохраняем результат, удаляем PDB со структурой и файл stride
- Путь для хранения вырезанных файлов также задаем в конфигурации. Они нужны только для расчета и не нужны впоследствии
- Учитываем, что целевой файл отсортирован по experiment_id (т.е. PDB ID). таким образом, читаем файл по несколько строк за раз, по PDB_ID; читаем один pdb-файл за раз, и держа этот файл в памяти, режем на фрагменты. Все мотивы в одном experiment_id обрабатываем в одном потоке
- Преобразование огромного числа получившихся JSON'ов в итоговые CSV-файлы делаем отдельным скриптом, который запускается отдельно - со своими параметрами START_INDEX, STOP_INDEX. В имени файла отражаем эти величины.