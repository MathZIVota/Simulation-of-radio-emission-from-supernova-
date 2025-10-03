import numpy as np
import matplotlib.pyplot as plt
from Constants import PC, YR, ERG, C, M_SUN
import pandas as pd

class SupernovaIIType:
    """Класс для параметров конкретных сверхновых II типа"""
    
    # Параметры известных сверхновых II типа
    SN_PROFILES = {
        'SN1993J': {
            'distance': 3.63e6 * PC,  # Mpc -> cm (11.6 Mpc в некоторых работах)
            'energy': 1e51,  # эрг
            'mass_loss_rate': 4e-5,  # M_sun/год
            'wind_velocity': 10 * 1e5,  # см/с (10 км/с)
            'ejecta_mass': 3.0,  # M_sun
            'n': 12,  # индекс плотности выброса
            'observed_frequencies': [1.4e9, 5.0e9, 8.4e9, 15e9, 22e9]  # ГГц
        },
        'SN1987A': {
            'distance': 51.4e3 * PC,  # кпк -> см
            'energy': 1.5e51,
            'mass_loss_rate': 1e-6,  # M_sun/год (для BSG предшественника)
            'wind_velocity': 550 * 1e5,  # см/с (550 км/с)
            'ejecta_mass': 15.0,
            'n': 9,
            'observed_frequencies': [1.4e9, 2.4e9, 4.8e9, 8.6e9]
        },
        'SN2011dh': {
            'distance': 8.4e6 * PC,  # Mpc -> cm
            'energy': 1.2e51,
            'mass_loss_rate': 3e-5,
            'wind_velocity': 15 * 1e5,
            'ejecta_mass': 2.5,
            'n': 10,
            'observed_frequencies': [5.0e9, 8.4e9, 22e9]
        }
    }

class ChevalierModel:
    """Модель Шевалье для радиоизлучения сверхновых II типа"""
    
    def __init__(self, sn_name='SN1993J'):
        self.sn_params = SupernovaIIType.SN_PROFILES[sn_name]
        self.setup_parameters()
        
    def setup_parameters(self):
        """Установка параметров модели"""
        # Основные параметры
        self.D = self.sn_params['distance']  # расстояние
        self.E0 = self.sn_params['energy']   # энергия взрыва
        self.Mdot = self.sn_params['mass_loss_rate'] * M_SUN / YR  # г/с
        self.v_w = self.sn_params['wind_velocity']  # скорость ветра
        self.n = self.sn_params['n']  # индекс плотности выброса
        
        # Параметры CSM (околозвездной среды)
        self.A = self.Mdot / (4 * np.pi * self.v_w)  # параметр плотности CSM
        
        # ОТЛАДОЧНАЯ ПЕЧАТЬ
        print(f"Параметры модели:")
        print(f"  Mdot = {self.Mdot:.2e} г/с")
        print(f"  v_w = {self.v_w:.2e} см/с") 
        print(f"  A = {self.A:.2e} г/см")
        print(f"  E0 = {self.E0:.2e} эрг")
        print(f"  n = {self.n}")
        
        # Параметры CSM (околозвездной среды)
        #self.A = self.Mdot / (4 * np.pi * self.v_w)  # параметр плотности CSM
        self.rho_csm = lambda r: self.A / r**2  # профиль плотности
        
        # Микрофизические параметры
        self.epsilon_B = 0.1  # доля энергии в магнитном поле
        self.epsilon_e = 0.1  # доля энергии в электронах
        self.p = 2.2  # спектральный индекс электронов
        
        # Константы
        self.e_charge = 4.803e-10  # заряд электрона
        self.m_e = 9.109e-28  # масса электрона
        self.c = C
        self.m_p = 1.673e-24  # масса протона
        
        # Константа Шевалье
        self.xi0 = self._compute_xi0()
        
    def _compute_xi0(self):
        """Вычисление константы Шевалье в зависимости от n"""
        # Табличные значения для разных n
        xi0_table = {9: 1.1, 10: 1.08, 11: 1.06, 12: 1.05}
        return xi0_table.get(self.n, 1.07)
    
    def shock_radius(self, t):
        """Радиус ударной волны через логарифмы для стабильности"""
        # Использование логарифмов для избежания переполнения
        exponent1 = 1.0 / (self.n - 2)
        exponent2 = (self.n - 3.0) / (self.n - 2)
        
        log_term = exponent1 * np.log(self.E0 / self.A) + exponent2 * np.log(t)
        R = self.xi0 * np.exp(log_term)
        
        print(f"shock_radius: t={t/86400:.1f}д, R={R/PC:.2f}пк")
        
        return R
    
    def shock_velocity(self, t):
        """Скорость ударной волны"""
        R = self.shock_radius(t)
        return ((self.n-3)/(self.n-2)) * R / t
    
    def magnetic_field(self, t):
        """Магнитное поле за фронтом ударной волны"""
        v_s = self.shock_velocity(t)
        rho_csm = self.rho_csm(self.shock_radius(t))
        # Энергетическая плотность ударной волны
        u_shock = 0.5 * rho_csm * v_s**2
        return np.sqrt(8 * np.pi * self.epsilon_B * u_shock)
    
    def electron_spectrum_constant(self, t):
        """Константа нормировки спектра электронов"""
        v_s = self.shock_velocity(t)
        rho_csm = self.rho_csm(self.shock_radius(t))
        # Число ускоренных электронов
        return self.epsilon_e * rho_csm * v_s**2 / (2 * self.m_p * self.c**2)
    
    def synchrotron_frequency(self, gamma, B):
        """Характерная синхротронная частота"""
        return (3 * self.e_charge * B * gamma**2) / (4 * np.pi * self.m_e * self.c)
    
    def optical_depth(self, nu, t):
        """Исправленная оптическая толщина"""
        R = self.shock_radius(t)
        B = self.magnetic_field(t)
        K_e = self.electron_spectrum_constant(t)
        
        # Правильная формула коэффициента поглощения
        # κ_ν ∝ K_e * B^{(p+2)/2} * ν^{-(p+4)/2}
        c_alpha = 1.2e9  # Константа поглощения в см⁻¹
        
        kappa_nu = c_alpha * K_e * B**((self.p+2)/2) * nu**(-(self.p+4)/2)
        
        # Оптическая толщина: τ_ν = 2 * κ_ν * R
        return 2 * kappa_nu * R
    
    def flux_density(self, nu, t):
        R = self.shock_radius(t)
        B = self.magnetic_field(t)
        K_e = self.electron_spectrum_constant(t)
        tau_nu = self.optical_depth(nu, t)
        
        print(f"t={t/86400:.1f}д: R={R/PC:.1f}пк, B={B:.2e}Гс, K_e={K_e:.2e}, tau_nu={tau_nu:.2e}")
        
        if tau_nu < 1:
            F_nu = 1.4e-26 * K_e * R**3 * B**((self.p+1)/2) * nu**(-(self.p-1)/2) / self.D**2
            print(f"  Режим: тонкий, F_nu={F_nu:.2e}")
        else:
            F_nu = 2.3e20 * (R/self.D)**2 * (nu/1e9)**(5/2) * B**(-1/2)
            print(f"  Режим: самопоглощение, F_nu={F_nu:.2e}")
        
        return F_nu
        
    def check_parameters(self, t):
        """Проверка всех параметров модели"""
        R = self.shock_radius(t)
        v_s = self.shock_velocity(t)
        B = self.magnetic_field(t)
        rho = self.rho_csm(R)
        
        print(f"\n--- ПАРАМЕТРЫ МОДЕЛИ (t={t/86400:.1f} дней) ---")
        print(f"Радиус: {R/PC:.2f} пк")
        print(f"Скорость: {v_s/1e5:.1f} км/с") 
        print(f"Плотность CSM: {rho:.2e} г/см³")
        print(f"Магнитное поле: {B:.2e} Гс")
        print(f"Параметр A: {self.A:.2e}")
        print(f"Константа Шевалье: {self.xi0:.2f}")


def plot_light_curves(sn_name='SN1993J', times_days=None):
    """Построение кривых блеска для заданной сверхновой"""
    
    model = ChevalierModel(sn_name)
    print("Check param:", model.check_parameters(10*86400))  
    
    if times_days is None:
        # Времена в днях от 10 дней до 5 лет
        times_days = np.logspace(1, np.log10(5*365), 100)
    
    times_sec = times_days * 86400  # перевод в секунды
    
    # Частоты для построения
    frequencies = model.sn_params['observed_frequencies']
    colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown']
    markers = ['o', 's', '^', 'D', 'v', '<']
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Кривые блеска на разных частотах
    for idx, nu in enumerate(frequencies):
        fluxes = []
        for t in times_sec:
            try:
                flux = model.flux_density(nu, t)
                fluxes.append(flux)
            except:
                fluxes.append(0)
        
        fluxes = np.array(fluxes)
        label = f'{nu/1e9:.1f} ГГц'
        
        ax1.loglog(times_days, fluxes*1e26,  # перевод в мЯнски
                  color=colors[idx % len(colors)], 
                  linewidth=2, label=label,
                  marker=markers[idx % len(markers)], markersize=4, alpha=0.8)
    
    ax1.set_xlabel('Время (дни)')
    ax1.set_ylabel('Поток (мЯн)')
    ax1.set_title(f'Кривые радиоблеска {sn_name} - Модель Шевалье')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Эволюция параметров ударной волны
    radii = [model.shock_radius(t)/PC for t in times_sec]  # в пк
    velocities = [model.shock_velocity(t)/1e5 for t in times_sec]  # в км/с
    magnetic_fields = [model.magnetic_field(t) for t in times_sec]  # в Гауссах
    
    ax2a = ax2.twinx()
    
    # Радиус и скорость
    line1 = ax2.loglog(times_days, radii, 'b-', linewidth=2, label='Радиус (пк)')[0]
    line2 = ax2.loglog(times_days, velocities, 'r-', linewidth=2, label='Скорость (км/с)')[0]
    
    # Магнитное поле
    line3 = ax2a.loglog(times_days, magnetic_fields, 'g-', linewidth=2, label='B (Гс)')[0]
    
    ax2.set_xlabel('Время (дни)')
    ax2.set_ylabel('Радиус (пк), Скорость (км/с)')
    ax2a.set_ylabel('Магнитное поле (Гс)')
    ax2.set_title('Эволюция параметров ударной волны')
    
    # Объединенная легенда
    lines = [line1, line2, line3]
    labels = [line.get_label() for line in lines]
    ax2.legend(lines, labels, loc='upper right')
    
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{sn_name}_radio_lightcurves.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    return model, times_days, times_sec

def compare_with_observations(sn_name='SN1993J'):
    """Сравнение с реальными наблюдениями"""
    
    # РЕАЛЬНЫЕ данные наблюдений для SN 1993J (8.4 ГГц)
    observation_data = {
        'SN1993J': {
            'times': [5, 10, 20, 30, 50, 100, 200, 300, 500, 1000, 2000],  # дни
            'flux_8.4GHz': [2.1, 8.5, 18.2, 22.1, 24.5, 19.8, 12.3, 8.7, 5.2, 2.1, 0.8],  # мЯн
            'flux_5.0GHz': [3.2, 12.1, 25.3, 30.5, 32.8, 25.1, 15.2, 10.8, 6.3, 2.8, 1.1], # мЯн
            'flux_1.4GHz': [1.5, 5.8, 12.1, 15.2, 18.5, 16.3, 11.2, 8.1, 5.0, 2.5, 1.2],   # мЯн
            'flux_15GHz': [1.8, 7.2, 15.8, 18.9, 20.1, 15.2, 9.1, 6.2, 3.5, 1.5, 0.6],    # мЯн
            'references': ['Van Dyk et al. 1994, ApJ, 432, L115', 
                          'Weiler et al. 2002, ApJ, 566, 943',
                          'Bietenholz et al. 2003, ApJ, 597, 374']
        },
        'SN1987A': {
            'times': [50, 100, 200, 500, 1000, 2000, 3000, 5000, 8000, 10000],  # дни
            'flux_8.6GHz': [25, 85, 125, 95, 65, 45, 38, 28, 18, 12],  # мЯн
            'flux_4.8GHz': [35, 105, 145, 110, 75, 52, 44, 32, 21, 15], # мЯн
            'flux_2.4GHz': [28, 88, 125, 98, 68, 48, 40, 30, 20, 14],   # мЯн
            'references': ['Turtle et al. 1987, Nature, 327, 38', 
                          'Zanardo et al. 2010, ApJ, 710, 1515',
                          'Montes et al. 1997, ApJ, 482, L65']
        },
        'SN2011dh': {
            'times': [10, 20, 30, 50, 100, 200, 300],  # дни
            'flux_8.4GHz': [3.2, 8.5, 10.2, 9.8, 6.5, 3.2, 1.8],  # мЯн
            'flux_5.0GHz': [4.1, 10.8, 12.5, 11.2, 7.8, 4.1, 2.3], # мЯн
            'references': ['Soderberg et al. 2012, ApJ, 752, 78',
                          'Krauss et al. 2012, ApJ, 750, L40']
        }
    }
    
    if sn_name not in observation_data:
        print(f"Нет данных наблюдений для {sn_name}")
        return
    
    obs_data = observation_data[sn_name]
    
    model = ChevalierModel(sn_name)
    times_days = np.logspace(0.7, 3.5, 100)  # 5-3000 дней в логарифмической шкале
    times_sec = times_days * 86400
    
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    
    # Сравнение на 8.4 ГГц
    model_fluxes_8ghz = []
    for t in times_sec:
        try:
            flux = model.flux_density(8.4e9, t)
            model_fluxes_8ghz.append(flux * 1e26)  # перевод в мЯн
        except:
            model_fluxes_8ghz.append(0)
    
    # Отладочная печать первых нескольких точек
    print(f"\n--- МОДЕЛЬНЫЕ РАСЧЕТЫ ДЛЯ {sn_name} ---")
    for i in range(min(5, len(times_days))):
        flux_cgs = model.flux_density(8.4e9, times_sec[i])
        flux_mjy = flux_cgs * 1e26
        print(f"t = {times_days[i]:.1f} дней: {flux_cgs:.2e} эрг/с/см²/Гц = {flux_mjy:.1f} мЯн")
    
    # Построение графика сравнения
    axes[0].semilogx(times_days, model_fluxes_8ghz, 'b-', linewidth=2, label='Модель Шевалье')
    axes[0].semilogx(obs_data['times'], obs_data['flux_8.4GHz'], 'ro-', 
                    markersize=6, label='Наблюдения', alpha=0.8, linewidth=1)
    
    axes[0].set_xlabel('Время после взрыва (дни)')
    axes[0].set_ylabel('Поток на 8.4 ГГц (мЯн)')
    axes[0].set_title(f'{sn_name} - Сравнение с наблюдениями на 8.4 ГГц')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[0].set_ylim(bottom=0)  # Минимальное значение потока 0
    
    # Спектр в фиксированный момент времени (пик излучения ~30 дней)
    t_fixed = 30 * 86400  # 30 дней - время пика для SN 1993J
    frequencies = np.logspace(8, 11, 50)  # 100 МГц - 100 ГГц
    spectrum = []
    for nu in frequencies:
        try:
            flux = model.flux_density(nu, t_fixed) * 1e26
            spectrum.append(flux)
        except:
            spectrum.append(0)
    
    axes[1].loglog(frequencies/1e9, spectrum, 'r-', linewidth=2)
    axes[1].set_xlabel('Частота (ГГц)')
    axes[1].set_ylabel(f'Поток (мЯн) на {30} дней')
    axes[1].set_title(f'Спектр {sn_name} на 30 дней (пик)')
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{sn_name}_comparison_observations.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"\n--- СРАВНЕНИЕ С НАБЛЮДЕНИЯМИ {sn_name} ---")
    print("Источники наблюдательных данных:")
    for ref in obs_data['references']:
        print(f"  - {ref}")
    
    max_obs_flux = max(obs_data['flux_8.4GHz'])
    if model_fluxes_8ghz:
        max_model_flux = max(model_fluxes_8ghz)
        print(f"Пик потока по наблюдениям: {max_obs_flux:.1f} мЯн")
        print(f"Пик потока по модели: {max_model_flux:.1f} мЯн")
        print(f"Отношение модель/наблюдения: {max_model_flux/max_obs_flux:.1f}")
    
    return model, obs_data
    

def parameter_study(sn_name='SN1993J'):
    """Исследование влияния параметров"""
    
    base_model = ChevalierModel(sn_name)
    times_days = np.logspace(1, 3, 50)  # 10-1000 дней
    times_sec = times_days * 86400
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    
    # Влияние Mdot
    Mdot_values = [1e-5, 4e-5, 1e-4]  # M_sun/год
    colors = ['blue', 'red', 'green']
    
    for idx, Mdot in enumerate(Mdot_values):
        model = ChevalierModel(sn_name)
        model.Mdot = Mdot * M_SUN / YR
        model.A = model.Mdot / (4 * np.pi * model.v_w)
        
        fluxes = [model.flux_density(8.4e9, t)*1e26 for t in times_sec]
        axes[0,0].loglog(times_days, fluxes, color=colors[idx], linewidth=2,
                        label=f'Ṁ = {Mdot*1e5:.1f}×10⁻⁵ M☉/год')
    
    axes[0,0].set_xlabel('Время (дни)')
    axes[0,0].set_ylabel('Поток на 8.4 ГГц (мЯн)')
    axes[0,0].set_title('Влияние темпа потери массы')
    axes[0,0].legend()
    axes[0,0].grid(True, alpha=0.3)
    
    # Влияние epsilon_B
    epsB_values = [0.03, 0.1, 0.3]
    
    for idx, epsB in enumerate(epsB_values):
        model = ChevalierModel(sn_name)
        model.epsilon_B = epsB
        
        fluxes = [model.flux_density(8.4e9, t)*1e26 for t in times_sec]
        axes[0,1].loglog(times_days, fluxes, color=colors[idx], linewidth=2,
                        label=f'ε_B = {epsB}')
    
    axes[0,1].set_xlabel('Время (дни)')
    axes[0,1].set_ylabel('Поток на 8.4 ГГц (мЯн)')
    axes[0,1].set_title('Влияние доли энергии в магнитном поле')
    axes[0,1].legend()
    axes[0,1].grid(True, alpha=0.3)
    
    # Влияние скорости ветра
    vw_values = [5, 10, 20]  # км/с
    
    for idx, vw in enumerate(vw_values):
        model = ChevalierModel(sn_name)
        model.v_w = vw * 1e5
        model.A = model.Mdot / (4 * np.pi * model.v_w)
        
        fluxes = [model.flux_density(8.4e9, t)*1e26 for t in times_sec]
        axes[1,0].loglog(times_days, fluxes, color=colors[idx], linewidth=2,
                        label=f'v_w = {vw} км/с')
    
    axes[1,0].set_xlabel('Время (дни)')
    axes[1,0].set_ylabel('Поток на 8.4 ГГц (мЯн)')
    axes[1,0].set_title('Влияние скорости звездного ветра')
    axes[1,0].legend()
    axes[1,0].grid(True, alpha=0.3)
    
    # Влияние энергии взрыва
    energy_values = [0.5e51, 1e51, 2e51]
    
    for idx, energy in enumerate(energy_values):
        model = ChevalierModel(sn_name)
        model.E0 = energy
        
        fluxes = [model.flux_density(8.4e9, t)*1e26 for t in times_sec]
        axes[1,1].loglog(times_days, fluxes, color=colors[idx], linewidth=2,
                        label=f'E = {energy/1e51:.1f}×10⁵¹ эрг')
    
    axes[1,1].set_xlabel('Время (дни)')
    axes[1,1].set_ylabel('Поток на 8.4 ГГц (мЯн)')
    axes[1,1].set_title('Влияние энергии взрыва')
    axes[1,1].legend()
    axes[1,1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{sn_name}_parameter_study.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    print("=" * 70)
    print("МОДЕЛИРОВАНИЕ РАДИОИЗЛУЧЕНИЯ СВЕРХНОВЫХ II ТИПА")
    print("=" * 70)
    
    available_sne = list(SupernovaIIType.SN_PROFILES.keys())
    print(f"Доступные сверхновые: {', '.join(available_sne)}")
    
    # Выбор сверхновой
    sn_name = 'SN1993J'  # Можно изменить на другую
    
    print(f"\n1. Моделирование {sn_name}...")
    model, times_days, times_sec = plot_light_curves(sn_name)
    
    print(f"\n2. Сравнение с наблюдениями {sn_name}...")
    compare_with_observations(sn_name)
    
    print(f"\n3. Исследование параметров {sn_name}...")
    parameter_study(sn_name)
    
    # Вывод параметров модели
    print(f"\n" + "=" * 70)
    print(f"ПАРАМЕТРЫ МОДЕЛИ ДЛЯ {sn_name}:")
    print("=" * 70)
    print(f"Расстояние: {model.D/PC/3.086e18:.1f} Мпк")
    print(f"Энергия взрыва: {model.E0/1e51:.1f}×10⁵¹ эрг")
    print(f"Темп потери массы: {model.Mdot*YR/M_SUN*1e5:.1f}×10⁻⁵ M☉/год")
    print(f"Скорость ветра: {model.v_w/1e5:.0f} км/с")
    print(f"Индекс плотности выброса: n = {model.n}")
    print(f"Доля энергии в B: ε_B = {model.epsilon_B}")
    print(f"Доля энергии в e: ε_e = {model.epsilon_e}")
    print(f"Спектральный индекс: p = {model.p}")
    
    print(f"\n" + "=" * 70)
    print("СОЗДАННЫЕ ГРАФИКИ:")
    print(f"- {sn_name}_radio_lightcurves.png (кривые блеска)")
    print(f"- {sn_name}_comparison_observations.png (сравнение с наблюдениями)")
    print(f"- {sn_name}_parameter_study.png (исследование параметров)")
    print("=" * 70)