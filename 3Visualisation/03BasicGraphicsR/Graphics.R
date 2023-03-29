# Визуализация данных стандартными средствами R ---------------------------

# Точечная диаграмма


# [, 1]	mpg	Miles/(US) gallon
# [, 2]	cyl	Number of cylinders
# [, 3]	disp	Displacement (cu.in.)
# [, 4]	hp	Gross horsepower
# [, 5]	drat	Rear axle ratio
# [, 6]	wt	Weight (1000 lbs)
# [, 7]	qsec	1/4 mile time
# [, 8]	vs	Engine (0 = V-shaped, 1 = straight)
# [, 9]	am	Transmission (0 = automatic, 1 = manual)
# [,10]	gear	Number of forward gears
# [,11]	carb	Number of carburetors

data_cars <- mtcars
str(data_cars)

data_cars <- transform(
    `_data` = data_cars,
    vs = as.factor(vs),
    gear = as.factor(gear),
    am = as.factor(am),
    carb = as.factor(carb),
    cyl = as.factor(cyl)
)

str(data_cars)

plot(data_cars, pch = 20)

attach(data_cars)
plot(
    x = drat, 
    y = mpg,
    pch = 20,
    col = gear,
    cex = I(1.5),
    main = "Количество передач переднего хода",
    xlab = "Передаточное отношение заднего моста",
    ylab = "Мили на галлон"
)
legend(
    x = min(drat),
    y = max(mpg),
    col = unique(gear),
    legend = unique(gear),
    pch = 20, 
    title = "Передач"
)
grid(
    
)
detach(data_cars)



# Коробки с усами ---------------------------------------------------------

# Длина зубов морских свинок
?ToothGrowth
data_to_plot <- ToothGrowth


attach(data_to_plot)
boxplot(
    formula = len ~ supp, 
    horizontal = FALSE, 
    notch = TRUE,
    data = data_to_plot,
    main = "Зависимость от вида витамина",
    xlab = "Вид витамина",
    ylab = "Длина зубов"
)
stripchart(
    len ~ supp,
    add = TRUE, 
    method = "jitter",
    pch = 19,
    col = 1:length(unique(supp)),
    vertical = TRUE
)
detach(data_to_plot)


attach(data_to_plot)
boxplot(
    formula = len ~ dose, 
    horizontal = FALSE, 
    notch = TRUE,
    data = data_to_plot,
    main = "Зависимость от дозы витамина",
    xlab = "Доза витамина, мг/день",
    ylab = "Длина зубов"
)
stripchart(
    len ~ dose,
    add = TRUE, 
    method = "jitter",
    pch = 19,
    col = 1:length(unique(dose)),
    vertical = TRUE
)
grid()
detach(data_to_plot)


# Столбчатые диаграммы ----------------------------------------------------

data_cars

data_agg <- aggregate(
    mpg ~ cyl + vs, 
    data = data_cars, 
    FUN = mean
)

attach(data_agg)
barplot(
    height = mpg, 
    names.arg = paste0("cyl:",cyl,", vs:", vs), 
    horiz = FALSE, 
    col = 1:nrow(data_agg), 
    main = "Средний расход топлива авто",
    sub = "В зависимости от расположения и количества цилиндров", 
    ylab = "Расход топлива, миль/галон" 
)
detach(data_agg)


# Гистограмма -------------------------------------------------------------
data_url <- "https://raw.githubusercontent.com/qwerty29544/RpracticeBook/master/2Data/01FlatTables/weatherAUS.csv"
weather_AUS <- read.csv(
    file = data_url, 
    sep = ",", 
    header = TRUE, 
    dec = ".", 
    encoding = "UTF-8"
)


Albury <- subset(weather_AUS, Location == "Albury")

attach(Albury)
x <- seq(
    from = min(Albury$MaxTemp, na.rm = T), 
    to = max(Albury$MaxTemp, na.rm = T), 
    length.out = 10000
)
f <- dnorm(
    x = x, 
    mean = mean(Albury$MaxTemp, na.rm = T),
    sd = sd(Albury$MaxTemp, na.rm = T)
)
hist(
    x = MaxTemp, 
    freq = F,
    col = "Azure"
)
lines(
    x = x, 
    y = f,
    col = "blue"
)
lines(
    density(MaxTemp, na.rm = T), 
    lwd = 2, 
    col = 'red'
)
detach(Albury)


