<script>
  import { language } from '#lib/stores/language.js';
  const options = ['en', 'fr'];
</script>

<div class="language-toggle" aria-label="Language" data-testid="language-toggle">
  {#each options as option}
    <button type="button" class:active={$language === option} aria-pressed={$language === option} onclick={() => language.set(option)}>
      {option.toUpperCase()}
    </button>
  {/each}
</div>

<style>
  .language-toggle{display:inline-flex;align-items:center;border:1px solid var(--border-subtle);border-radius:6px;overflow:hidden;background:var(--bg-tertiary)}
  button{border:0;background:transparent;color:var(--text-secondary);cursor:pointer;font:inherit;font-family:var(--font-mono);font-size:var(--text-xs);font-weight:700;padding:var(--space-2) var(--space-3)}
  button:hover{color:var(--text-primary);background:var(--bg-secondary)}
  button.active{background:var(--text-primary);color:var(--bg-primary)}
</style>
