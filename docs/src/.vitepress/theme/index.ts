// The DocumenterVitepress default theme, which docs/make.jl copies to
// base.ts from the installed package so it tracks that package's version,
// extended with the large-table enhancement.
import { onMounted } from 'vue'
import type { Theme as ThemeConfig } from 'vitepress'
import BaseTheme from './base'
import { watchTables } from './tables'
import './tables.css'
// __DV_PLUGIN_COMPONENT_IMPORTS__

export const Theme: ThemeConfig = {
  ...BaseTheme,
  enhanceApp(ctx) {
    BaseTheme.enhanceApp?.(ctx)
    const { app } = ctx
    // __DV_PLUGIN_COMPONENT_REGISTRATIONS__
  },
  setup() {
    BaseTheme.setup?.()
    // After hydration, so the server-rendered tables match on the client.
    onMounted(watchTables)
  },
}
export default Theme
