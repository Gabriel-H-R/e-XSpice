/**
 * @file e_XSpice_rcp.h
 * @brief Implementação de Ponteiros Inteligentes de Alta Performance para o CAS.
 * * Este arquivo define a fundação de gerenciamento de memória baseada no design
 * da SymEngine. Ele provê ponteiros intrusivos para máxima performance (Cache locality
 * e alocação única) e um modo de Debug para rastreamento de vazamentos de memória.
 */

#pragma once

#include <atomic>
#include <utility>
#include <cassert>
#include <type_traits>
#include <stdexcept>


// Suporte retroativo para C++20 [[unlikely]] / [[likely]]
#if defined(__has_cpp_attribute) && __has_cpp_attribute(unlikely)
    #define SYMENGINE_UNLIKELY [[unlikely]]
#else
    #define SYMENGINE_UNLIKELY
#endif


// ============================================================================
// CONFIGURAÇÃO DO MODO DE DEBUG
// ============================================================================
// Na SymEngine original, o modo de debug usa a biblioteca Teuchos. 
// Na nossa implementação, podemos ativar verificações extras pesadas ou 
// cair para o std::shared_ptr definindo a macro abaixo durante a compilação.
// #define SYMENGINE_DEBUG_MODE 

#ifdef SYMENGINE_DEBUG_MODE
    // Macros de asserção para facilitar a detecção de ponteiros nulos
    #define SYMENGINE_ASSERT(cond, msg) \
        if (!(cond)) { throw std::runtime_error(msg); }
#else
    // Em modo release, as asserções viram código nulo (custo zero)
    #define SYMENGINE_ASSERT(cond, msg) assert(cond)
#endif


/**
 * @brief Classe base para habilitação de contagem de referências intrusiva.
 * * Qualquer classe do sistema de álgebra (ex: Symbol, Add, Mul) deve herdar
 * desta classe. O contador de referência é alocado adjacente aos dados do objeto
 * garantindo localidade de cache e velocidade extrema na manipulação do DAG matemático.
 */
class RefCounted {
private:
    /// @brief Contador atômico seguro para paralelismo (Thread-safe).
    mutable std::atomic<unsigned int> refcount_;

public:
    /**
     * @brief Construtor padrão. Inicializa a contagem em 0.
     */
    constexpr RefCounted() noexcept : refcount_(0) {}
    
    /**
     * @brief Destrutor virtual padrão.
     * @note Essencial para garantir a deleção polimórfica correta das classes derivadas.
     */
    virtual ~RefCounted() = default;

    /**
     * @brief Construtor de cópia anulado para referências.
     * @note Se um objeto é copiado, o novo objeto começa com 0 referências.
     */
    RefCounted(const RefCounted&) noexcept : refcount_(0) {}
    RefCounted& operator=(const RefCounted&) noexcept { return *this; }

    /**
     * @brief Incrementa a contagem de referência.
     * @note Utiliza memory_order_relaxed para evitar barreiras de sincronização na CPU,
     * maximizando a velocidade já que a ordem de incremento não afeta outras threads.
     */
    inline void add_ref() const noexcept {
        refcount_.fetch_add(1, std::memory_order_relaxed);
    }

    /**
     * @brief Decrementa a contagem de referência e deleta o objeto se chegar a zero.
     * @note Utiliza memory_order_acq_rel para garantir que todas as escritas no objeto
     * por qualquer thread sejam finalizadas antes da destruição (delete this).
     */
    inline void dec_ref() const noexcept {
        // memory_order_release garante que modificações neste objeto
        // feitas por esta thread sejam concluídas antes do decremento.
        if (refcount_.fetch_sub(1, std::memory_order_release) == 1) SYMENGINE_UNLIKELY {
            // Apenas a ultima thread (que deletará o objeto) paga o preço do acquire.
            // Sincroniza com as operações 'release' de todas as outras threads.
            std::atomic_thread_fence(std::memory_order_acquire);
            delete this;
        }
    }
    // inline void dec_ref() const noexcept {
    //     if (refcount_.fetch_sub(1, std::memory_order_acq_rel) == 1) {
    //         delete this;
    //     }
    // }


    /**
     * @brief Retorna a quantidade de referências ativas apontando para este objeto.
     * @return unsigned int Número de referências.
     */
    [[nodiscard]] inline unsigned int get_refcount() const noexcept {
        return refcount_.load(std::memory_order_relaxed);
    }
};


/**
 * @brief Ponteiro inteligente de contagem de referência (Reference Counted Pointer).
 * * @tparam T O tipo da classe armazenada (deve herdar de RefCounted).
 * * Gerencia automaticamente o ciclo de vida dos nós matemáticos. Implementa 
 * 'Move Semantics' (C++11) para transferências de propriedade a custo zero.
 */
template <typename T>
class RCP {
private:
    T* ptr_; ///< Ponteiro bruto encapsulado.

    // Permite que RCP<Basic> leia o ponteiro privado de um RCP<Integer>
    template <typename U> friend class RCP;

    /**
     * @brief Decrementa a referência atual e zera o ponteiro.
     */
    inline void release() noexcept {
        if (ptr_ != nullptr) {
            ptr_->dec_ref();
            ptr_ = nullptr;
        }
    }

public:
    /** @brief Construtor nulo. */
    constexpr RCP() noexcept : ptr_(nullptr) {}
    constexpr RCP(std::nullptr_t) noexcept : ptr_(nullptr) {}

    /**
     * @brief Assume a posse de um ponteiro bruto recém-alocado.
     * @param p Ponteiro bruto para o objeto derivado de RefCounted.
     */
    explicit RCP(T* p) noexcept : ptr_(p) {
        if (ptr_) ptr_->add_ref();
    }

    /**
     * @brief Construtor de Cópia.
     * @param other Outro RCP a ser compartilhado.
     */
    RCP(const RCP& other) noexcept : ptr_(other.ptr_) {
        if (ptr_) ptr_->add_ref();
    }

    /**
     * @brief Construtor de Movimento.
     * @param other Outro RCP a ter sua posse "roubada" (custo zero na variável atômica).
     */
    RCP(RCP&& other) noexcept : ptr_(other.ptr_) {
        other.ptr_ = nullptr;
    }

    /** @brief Destrutor. Libera o nó do grafo. */
    ~RCP() { release(); }

    /**
     * @brief Operador de atribuição de cópia.
     * @param other O ponteiro a ser copiado.
     * @return RCP& Referência a si próprio.
     */
    // Copy Assignment otimizado
    RCP& operator=(const RCP& other) noexcept {
        if (ptr_ != other.ptr_) {
            // Incrementa o novo ANTES de decrementar o antigo (seguro contra auto-atribuição disfarçada)
            if (other.ptr_) other.ptr_->add_ref();
            T* old = ptr_;
            ptr_ = other.ptr_;
            // O desvio fica apenas aqui no final. Muito mais rápido para o pipeline.
            if (old) old->dec_ref();
        }
        return *this;
    }
    // RCP& operator=(const RCP& other) noexcept {
    //     if (this != &other) {
    //         if (ptr_ != other.ptr_) {
    //             release();
    //             ptr_ = other.ptr_;
    //             if (ptr_) ptr_->add_ref();
    //         }
    //     }
    //     return *this;
    // }

    /**
     * @brief Operador de atribuição de movimento.
     */
    //Move Assignment totalmente branchless no caminho principal
    RCP& operator=(RCP&& other) noexcept {
    if (this != &other) {
        T* old = ptr_;
        ptr_ = other.ptr_;
        other.ptr_ = nullptr;
        if (old) old->dec_ref();
    }
    return *this;
}
    // RCP& operator=(RCP&& other) noexcept {
    //     if (this != &other) {
    //         release();
    //         ptr_ = other.ptr_;
    //         other.ptr_ = nullptr;
    //     }
    //     return *this;
    // }

    /**
     * @brief Converte implicitamente cópias de RCP derivadas para RCP base.
     * @note Ex: RCP<Integer> para RCP<Basic>. Incrementa a referência.
     */
    template <typename U, typename std::enable_if<std::is_convertible<U*, T*>::value, int>::type = 0>
    RCP(const RCP<U>& other) noexcept : ptr_(other.ptr_) {
        if (ptr_) ptr_->add_ref();
    }

    /**
     * @brief Converte implicitamente movimentos de RCP derivadas para base.
     * @note Isso é exatamente o que é chamado no seu 'return integer(...)'.
     * Custo zero! O ponteiro é transferido sem tocar nas variáveis atômicas!
     */
    template <typename U, typename std::enable_if<std::is_convertible<U*, T*>::value, int>::type = 0>
    RCP(RCP<U>&& other) noexcept : ptr_(other.ptr_) {
        other.ptr_ = nullptr; // O temporário perde a posse
    }

    /**
     * @brief Retorna o ponteiro bruto subjacente.
     * @return T* O ponteiro bruto (pode ser nulo).
     */
    [[nodiscard]] inline T* get() const noexcept { return ptr_; }

    /**
     * @brief Desreferencia o ponteiro.
     * @note No modo de debug, checa agressivamente por nullptrs.
     * @return T& Referência ao objeto.
     */
    [[nodiscard]] inline T& operator*() const {
        SYMENGINE_ASSERT(ptr_ != nullptr, "Erro: Desreferenciamento de RCP nulo!");
        return *ptr_;
    }

    /**
     * @brief Acesso aos membros do ponteiro.
     * @return T* Ponteiro para o objeto.
     */
    [[nodiscard]] inline T* operator->() const {
        SYMENGINE_ASSERT(ptr_ != nullptr, "Erro: Acesso a membro em RCP nulo!");
        return ptr_;
    }

    /** @brief Avalia booleano (true se não for nulo). */
    explicit operator bool() const noexcept {
        return ptr_ != nullptr;
    }

    /** @brief Permite acesso rápido para trocar ponteiros. */
    void swap(RCP& other) noexcept { std::swap(ptr_, other.ptr_); }
};


/**
 * @brief Ponteiro não-dono leve para passagem de argumentos.
 * * @tparam T O tipo da classe apontada.
 * * Equivale ao ponteiro bruto tradicional C++, usado nas assinaturas
 * de funções internas onde não desejamos alterar o contador atômico do nó,
 * mas apenas inspecionar a expressão matemática (evitando gargalo de CPU).
 */
template <typename T>
class Ptr {
private:
    T* ptr_; ///< Ponteiro bruto encapsulado.

public:
    constexpr Ptr() noexcept : ptr_(nullptr) {}
    constexpr Ptr(std::nullptr_t) noexcept : ptr_(nullptr) {}
    
    /** @brief Construtor implícito a partir de um ponteiro bruto. */
    Ptr(T* p) noexcept : ptr_(p) {}

    /** * @brief Construtor que permite injetar um RCP diretamente como argumento Ptr. 
     * @param rcp O ponteiro inteligente dono.
     */
    Ptr(const RCP<T>& rcp) noexcept : ptr_(rcp.get()) {}

    [[nodiscard]] inline T* get() const noexcept { return ptr_; }
    [[nodiscard]] inline T& operator*() const {
        SYMENGINE_ASSERT(ptr_ != nullptr, "Erro: Desreferenciamento de Ptr nulo!");
        return *ptr_;
    }
    [[nodiscard]] inline T* operator->() const {
        SYMENGINE_ASSERT(ptr_ != nullptr, "Erro: Acesso em Ptr nulo!");
        return ptr_;
    }
    explicit operator bool() const noexcept { return ptr_ != nullptr; }
};

// ============================================================================
// UTILITÁRIOS GLOBAIS DE CRIAÇÃO E CASTING
// ============================================================================

/**
 * @brief Fábrica otimizada para alocação de classes matemáticas no RCP.
 * * Funciona como std::make_shared. Encaminha perfeitamente (Perfect Forwarding)
 * os argumentos para o construtor da classe T.
 * * @tparam T Tipo a ser criado.
 * @tparam Args Tipos dos argumentos do construtor.
 * @param args Argumentos a serem passados ao construtor de T.
 * @return RCP<T> Ponteiro RCP detentor do novo objeto.
 */
template <typename T, typename... Args>
[[nodiscard]] inline RCP<T> make_rcp(Args&&... args) {
    return RCP<T>(new T(std::forward<Args>(args)...));
}

/**
 * @brief Cast estático hiper-rápido para converter entre as hierarquias RCP.
 * @note Use apenas quando tiver certeza do tipo matemático alvo na hierarquia.
 * * @tparam T O tipo destino.
 * @tparam U O tipo origem.
 * @param rcp O ponteiro original.
 * @return RCP<T> O ponteiro convertido.
 */
template <typename T, typename U>
[[nodiscard]] inline RCP<T> rcp_static_cast(const RCP<U>& rcp) noexcept {
    return RCP<T>(static_cast<T*>(rcp.get()));
}

/**
 * @brief Cast dinâmico seguro em runtime.
 * @note Usa RTTI do C++. Se a conversão for inválida, retorna um RCP nulo.
 * * @tparam T O tipo destino.
 * @tparam U O tipo origem.
 * @param rcp O ponteiro original.
 * @return RCP<T> O ponteiro convertido ou nulo se falhar.
 */
template <typename T, typename U>
[[nodiscard]] inline RCP<T> rcp_dynamic_cast(const RCP<U>& rcp) noexcept {
    return RCP<T>(dynamic_cast<T*>(rcp.get()));
}