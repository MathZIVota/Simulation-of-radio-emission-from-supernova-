import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Константы
M_SUN = 1.988e33      # граммы
YR = 3.154e7          # секунды
PC = 3.085e18         # сантиметры (1 парсек)
E_CHARGE = 4.803e-10  # э.с.е. (заряд электрона)
M_E = 9.109e-28       # масса электрона, г
C = 3e10              # скорость света, см/с
M_P = 1.673e-24       # масса протона, г
SMALL_NUM = 1e-12     # Для обработки нулей


class ChevalierModel:
    """Модель Шевалье для радиоизлучения сверхновых II типа"""

    SN_PROFILES = {  # Определение SN_PROFILES как атрибута класса
        'SN1993J': {
            'distance': 3.63e24,       # cm (11.6 Mpc)
            'energy': 1e51,            # erg
            'mass_loss_rate': 4e-5,    # M☉/год
            'wind_velocity': 10e5,     # cm/с (10 km/s)
            'n': 10,  
            'initial_peak_flux': 20    # Приблизительная плотность потока в день пика (mJy)
        }
    }

    def __init__(self, sn_name='SN1993J'):
        self.sn_params = self.SN_PROFILES[sn_name]
        self.setup_parameters()

    def setup_parameters(self):
        """Установка параметров модели"""
        # Основные параметры
        self.D = self.sn_params['distance']  # расстояние в см
        self.E0 = self.sn_params['energy']   # энергия взрыва в эрг
        self.Mdot = self.sn_params['mass_loss_rate'] * M_SUN / YR  # г/с
        self.v_w = self.sn_params['wind_velocity']  # скорость ветра в см/с
        self.n = self.sn_params['n']  # индекс плотности выброса (n=10 для SN1993J)

        # Параметры CSM (околозвездной среды)
        self.A = self.Mdot / (4 * np.pi * self.v_w)  # параметр плотности CSM (g/cm)
        self.rho_csm = lambda r: self.A / r**2  # профиль плотности

        # Микрофизические параметры (типичные значения, их можно будет подгонять)
        self.epsilon_B = 0.1  # доля энергии в магнитном поле
        self.epsilon_e = 0.1  # доля энергии в электронах
        self.p = 2.2           # спектральный индекс электронов

        # Параметры для масштабирования
        self.t_peak_days = 100   # Приблизительное время пика в днях
        self.F_peak_MJY = self.sn_params['initial_peak_flux'] # Пиковый поток в mJy

        # Константа Шевалье
        self.xi0 = self._compute_xi0()

    def _compute_xi0(self):
        """Вычисление константы Шевалье в зависимости от n"""
        xi0_table = {9: 1.1, 10: 1.08, 11: 1.06, 12: 1.05}
        return xi0_table.get(self.n, 1.07)

    def shock_radius(self, t_sec):
        """Радиус ударной волны (t в секундах)"""
        # Формула для n=10
        # Используем логарифмическую форму для численной стабильности
        exponent1 = 1.0 / (self.n - 2)
        exponent2 = (self.n - 3.0) / (self.n - 2)

        # Избегаем log(0)
        if self.A <= 0 or self.E0 <= 0 or t_sec <= 0:
            return 0.0

        log_term = exponent1 * np.log(self.E0 / self.A) + exponent2 * np.log(t_sec)
        R = self.xi0 * np.exp(log_term)
        return R

    def shock_velocity(self, t_sec):
        """Скорость ударной волны (t в секундах)"""
        R = self.shock_radius(t_sec)
        if t_sec <= 0:
            return 0.0
        return ((self.n - 3) / (self.n - 2)) * R / t_sec

    def magnetic_field(self, t_sec):
        """Магнитное поле за фронтом ударной волны (приближение)"""
        R = self.shock_radius(t_sec)
        v_s = self.shock_velocity(t_sec)
        rho_csm = self.rho_csm(R)

        # Энергетическая плотность, которая сжимается
        u_shock = 0.5 * rho_csm * v_s**2

        # B^2 ~ epsilon_B * u_shock
        B = np.sqrt(8 * np.pi * self.epsilon_B * u_shock)
        return B

    def electron_spectrum_constant(self, t_sec):
        """Константа нормировки спектра электронов K_e"""
        v_s = self.shock_velocity(t_sec)
        rho_csm = self.rho_csm(self.shock_radius(t_sec))

        # K_e ~ epsilon_e * rho_csm * v_s^2 / m_p
        K_e = self.epsilon_e * rho_csm * v_s**2 / M_P
        return K_e

    def optical_depth(self, nu, t_sec):
        """Приближенная оптическая толщина"""
        R = self.shock_radius(t_sec)
        B = self.magnetic_field(t_sec)
        K_e = self.electron_spectrum_constant(t_sec)

        # Правильная формула коэффициента поглощения (Toro)
        c_alpha = 1.2e9 # Константа поглощения

        # kappa_nu ~ K_e * B^((p+2)/2) * nu^(-(p+4)/2)
        kappa_nu = c_alpha * K_e * B**((self.p+2)/2) * nu**(-(self.p+4)/2)

        tau_nu = 2 * kappa_nu * R
        return tau_nu

    def flux_density(self, nu, t_sec):
        """Вычисление плотности потока (в мЯн, 10^-29 эрг/с/см^2/Гц)"""
        R = self.shock_radius(t_sec)
        B = self.magnetic_field(t_sec)
        K_e = self.electron_spectrum_constant(t_sec)
        tau_nu = self.optical_depth(nu, t_sec)

        # 1 mJy = 1e-26 эрг/с/см^2/Гц в СГС единицах

        if tau_nu < 1:
            # Режим тонкого радиоисточника (self-absorption-free regime)
            # F_nu ~ K_e * R^3 * B^((p+1)/2) * nu^(-(p-1)/2) / D^2
            F_nu_cgs = 1.4e-26 * K_e * R**3 * B**((self.p+1)/2) * nu**(-(self.p-1)/2) / self.D**2

        else:
            # Режим самопоглощения (Self-Absorption regime)
            # F_nu ~ R^2 * nu^(5/2) * B^(-1/2)
            F_nu_cgs = 2.3e-20 * (R/self.D)**2 * (nu/1e9)**(5/2) * B**(-1/2)

        # Перевод из СГС (эрг/с/см^2/Гц) в мЯн (mJy)
        F_nu_mJy = F_nu_cgs * 1e3 # умножаем на 1e3, чтобы перевести в mJy
        return F_nu_mJy

    def check_parameters(self, t_sec):
        """Проверка всех параметров модели"""
        R = self.shock_radius(t_sec)
        v_s = self.shock_velocity(t_sec)
        B = self.magnetic_field(t_sec)
        rho = self.rho_csm(R)

        t_days = t_sec / 86400.0
        print(f"\n--- ПАРАМЕТРЫ МОДЕЛИ (t={t_days:.1f} дней) ---")
        print(f"Радиус (R): {R/PC:.2e} пк")
        print(f"Скорость (v_s): {v_s/1e5:.1f} км/с")
        print(f"Плотность CSM (rho): {rho:.2e} г/см^3")
        print(f"Магнитное поле (B): {B:.2e} Гс")
        print(f"Параметр A: {self.A:.2e}")
        print(f"Константа Шевалье (xi0): {self.xi0:.2f}")


# --- Наблюдательные данные SN 1993J на 8.4 GHz (Weiler et al. 2007) ---
#  ВРЕМЯ (дни) ПЛОТНОСТЬ ПОТОКА (мЯн) СТАНДАРТНАЯ ОШИБКА (мЯн)
OBS_TIME_DAYS = np.array([
    9.81, 14.82, 20.08, 34.74, 39.75, 53.83, 61.16, 74.74, 95.06, 115.25,
    135.50, 168.33, 199.83, 239.83, 269.83, 330.16, 369.83, 440.16, 539.83, 639.83, 740.16,
    840.16, 940.16
])

OBS_FLUX_MJY = np.array([
    2.05, 4.20, 6.20, 10.40, 11.60, 15.40, 16.30, 18.80, 20.80, 19.10,
    17.90, 15.30, 13.70, 11.60, 9.50, 7.70, 6.80, 5.20, 3.80, 3.00, 2.40,
    1.90, 1.50
])

OBS_ERROR_MJY = np.array([
    0.30, 0.30, 0.50, 0.70, 0.70, 0.90, 0.90, 1.10, 1.20, 1.10,
    1.10, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.20, 0.10,
    0.10, 0.10
])


# --- Определение функции, которая строит кривую Шевалье ---
def model_flux(t, t_peak=100, F_peak=20, alpha=0.8, beta=0.6):
    """Простая модель Шевалье (степенные законы)"""
    flux = np.zeros_like(t)
    # Рост до пика
    mask = t <= t_peak
    flux[mask] = F_peak * (t[mask]/t_peak)**alpha
    # Спад после пика
    flux[~mask] = F_peak * (t[~mask]/t_peak)**(-beta)
    return flux

# --- Определение функции, которая считает MSE для подгонки параметров ---
def calculate_mse(params, t_obs, F_obs, sigma_obs):
    """Вычисление среднеквадратичной ошибки MSE"""
    t_peak, F_peak, alpha, beta = params
    F_model = model_flux(t_obs, t_peak, F_peak, alpha, beta)
    mse = np.sum(((F_obs - F_model) / sigma_obs)**2) / len(t_obs)
    return mse


# --- ПОДГОНКА ПАРАМЕТРОВ ---
# Начальные приближения и границы для параметров
initial_guess = [100, 20, 0.8, 0.6]  # t_peak, F_peak, alpha, beta
bounds = (
    [50, 10, 0.0, 0.0], # Нижние границы
    [200, 30, 2.0, 2.0]  # Верхние границы
)


# --- Адаптация под curve_fit (для корректной обработки ошибок) ---
def model_curvefit(t, t_peak, F_peak, alpha, beta):
    return model_flux(t, t_peak, F_peak, alpha, beta)

# --- curve_fit ---
popt, pcov = curve_fit(
    model_curvefit,
    OBS_TIME_DAYS,
    OBS_FLUX_MJY,
    p0=initial_guess,
    sigma=OBS_ERROR_MJY,
    absolute_sigma=True,  # Учитываем, что sigma_obs задана абсолютно
    bounds=bounds
)


# Извлечение оптимальных параметров и их ошибок
t_peak_opt, F_peak_opt, alpha_opt, beta_opt = popt
t_peak_err, F_peak_err, alpha_err, beta_err = np.sqrt(np.diag(pcov)) # Квадратные корни диагональных элементов

# Создаем инстанс модели
model = ChevalierModel()

# Обновляем параметры модели (опционально)
model.sn_params['initial_peak_flux'] = F_peak_opt
model.t_peak_days = t_peak_opt
model.setup_parameters() # Пересчет параметров

# --- Строим график ---
t_model = np.linspace(5, 1200, 500)
F_model_opt = model_flux(t_model, t_peak_opt, F_peak_opt, alpha_opt, beta_opt)

# --- Выводим полученные параметры ---
print(f"Оптимальные параметры, полученные в результате подгонки:")
print(f"t_peak: {t_peak_opt:.2f} +/- {t_peak_err:.2f} дней")
print(f"F_peak: {F_peak_opt:.2f} +/- {F_peak_err:.2f} mJy")
print(f"alpha: {alpha_opt:.2f} +/- {alpha_err:.2f}")
print(f"beta: {beta_opt:.2f} +/- {beta_err:.2f}")

# Выводим финальный график
plt.figure(figsize=(10, 6))
plt.errorbar(OBS_TIME_DAYS, OBS_FLUX_MJY, yerr=OBS_ERROR_MJY, fmt='o', color='k',
            markersize=6, capsize=3, label='Weiler et al. (2007), 8.4 GHz')

plt.plot(t_model, F_model_opt, 'r-', lw=2, label='Подгонка (Chevalier Model)')
plt.axvline(t_peak_opt, color='gray', linestyle='--', alpha=0.7, label=f't_peak={t_peak_opt:.0f} дн.')
plt.xlabel('Время после взрыва (дни)', fontsize=12)
plt.ylabel('Плотность потока (мЯн)', fontsize=12)
plt.title('Кривая радиоизлучения SN1993J (8.4 GHz)', fontsize=14)
plt.yscale('log')
plt.xscale('log')
plt.legend(fontsize=10)
plt.grid(True, which='both', linestyle=':')
plt.tight_layout()
plt.show()